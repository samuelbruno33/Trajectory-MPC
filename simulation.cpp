#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "MPC.hpp"
#include "geometry.hpp"

using namespace std;

// --- UTILS ---
double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// --- PLANT (Modello Fisico Reale) ---
struct State { double x, y, theta, vx, vy, omega; };

State update_vehicle_physics(State s, double Fx, double delta, double dt) {
    double m = 320.0; double Iz = 98.03; double la = 0.792; double lb = 0.758;
    double B = 10.0, C = 1.9, D = 1.0; 

    // Limiti Fisici Reali
    if (Fx > 1500) Fx = 1500; if (Fx < -1500) Fx = -1500;
    double max_steer = 0.42; // 24 gradi
    if (delta > max_steer) delta = max_steer; if (delta < -max_steer) delta = -max_steer;

    double vx_safe = std::max(s.vx, 1.0);
    double alpha_f = delta - atan2(s.vy + s.omega * la, vx_safe);
    double alpha_r = -atan2(s.vy - s.omega * lb, vx_safe);

    double Fyf = D * sin(C * atan(B * alpha_f)) * m * 9.81 * (lb / (la + lb)); 
    double Fyr = D * sin(C * atan(B * alpha_r)) * m * 9.81 * (la / (la + lb));

    double dot_vx = (Fx - Fyf * sin(delta)) / m + s.vy * s.omega;
    double dot_vy = (Fyf * cos(delta) + Fyr) / m - s.vx * s.omega;
    double dot_omega = (la * Fyf * cos(delta) - lb * Fyr) / Iz;

    State next = s;
    next.vx += dot_vx * dt; next.vy += dot_vy * dt; next.omega += dot_omega * dt;
    next.x += (s.vx * cos(s.theta) - s.vy * sin(s.theta)) * dt;
    next.y += (s.vx * sin(s.theta) + s.vy * cos(s.theta)) * dt;
    next.theta = normalize_angle(s.theta + s.omega * dt);
    return next;
}

// --- TRACK GENERATION (Chicane: Dritto -> SX -> DX -> Dritto) ---
struct TrackPoint { double x, y, v_target; };

vector<TrackPoint> generate_track() {
    vector<TrackPoint> t;
    double ds = 0.5; // Risoluzione punti

    // 1. Primo Rettilineo (0 -> 80m) - ACCELERAZIONE
    for (double x = 0; x < 80.0; x += ds) {
        t.push_back({x, 0.0, 25.0}); // Target 90 km/h
    }

    // 2. Curva a Sinistra (80 -> 120m)
    // Usiamo un arco di cerchio o una coseno per fare la chicane
    // Spostamento laterale verso Y positivo (+4m)
    for (double x = 80.0; x < 120.0; x += ds) {
        double progress = (x - 80.0) / 40.0; // 0 to 1
        double y = 2.0 * (1.0 - cos(progress * M_PI)); // 0 -> 4m
        t.push_back({x, y, 12.0}); // Frenata a 43 km/h
    }

    // 3. Curva a Destra (120 -> 160m)
    // Ritorno a Y=0
    for (double x = 120.0; x < 160.0; x += ds) {
        double progress = (x - 120.0) / 40.0;
        double y = 4.0 - 2.0 * (1.0 - cos(progress * M_PI)); // 4m -> 0m
        t.push_back({x, y, 12.0}); // Mantieni bassa velocità
    }

    // 4. Secondo Rettilineo (160 -> 300m) - RI-ACCELERAZIONE
    for (double x = 160.0; x <= 300.0; x += ds) {
        t.push_back({x, 0.0, 25.0}); 
    }
    
    return t;
}

pair<vector<Point>, double> get_horizon(State s, const vector<TrackPoint>& track, int op) {
    double min_d = 1e9; int idx = 0;
    for(int i=0; i<track.size(); i++) {
        double d = sqrt(pow(s.x - track[i].x, 2) + pow(s.y - track[i].y, 2));
        if(d < min_d) { min_d = d; idx = i; }
    }
    
    vector<Point> wp;
    // Lookahead velocità intelligente: guarda 20 passi avanti per frenare in tempo
    int v_idx = std::min(idx + 30, (int)track.size()-1); 
    double v_ref = track[v_idx].v_target;

    for(int i=0; i<op; i++) {
        int ii = std::min(idx+i, (int)track.size()-1);
        wp.push_back(Point(track[ii].x, track[ii].y));
    }
    return {wp, v_ref};
}

// --- CONTROLLORI PID ---
class PID {
    double kp, ki, kd, I=0, prev=0, limit;
public:
    PID(double p, double i, double d, double lim) : kp(p), ki(i), kd(d), limit(lim) {}
    double calc(double ref, double meas, double dt) {
        double err = ref - meas;
        I += err*dt;
        if(I > 500) I=500; if(I<-500) I=-500; 
        double der = (err - prev)/dt;
        prev = err;
        double out = kp*err + ki*I + kd*der;
        return std::max(-limit, std::min(out, limit));
    }
    void reset() { I=0; prev=0; }
};

class PurePursuit {
public:
    double calc(State s, const vector<Point>& path) {
        // Lookahead meno aggressivo e più realistico
        double ld = std::clamp(0.5 * s.vx + 1.5, 2.5, 10.0);

        // Guadagno di sterzo ridotto (simula limiti fisici / understeer)
        const double wheelbase = 1.55;
        const double steer_gain = 0.85;

        for (size_t i = 0; i < path.size(); ++i) {
            double dx = path[i].x - s.x;
            double dy = path[i].y - s.y;
            double dist = std::hypot(dx, dy);

            if (dist > ld) {
                double alpha = atan2(dy, dx) - s.theta;
                alpha = normalize_angle(alpha);

                double delta = atan2(2.0 * wheelbase * sin(alpha), ld);
                return steer_gain * delta;
            }
        }
        return 0.0;
    }
};


// --- MAIN ---
int main() {
    // Check Critico
    ifstream check("cornering_stiffness_vs_vertical_load.txt");
    if(!check.good()) {
        cerr << "ERRORE: File 'cornering_stiffness_vs_vertical_load.txt' NON TROVATO!" << endl;
        cerr << "Assicurati che sia nella cartella build/Debug o build/" << endl;
        return -1;
    }

    auto track = generate_track();
    int op = 40; 
    double dt = 0.05;
    int steps = 500; // 25 secondi

    MPC mpc; 
    State s_mpc = {0,0,0, 1.0, 0,0};
    
    PID pid_gas(400, 5, 20, 1500); 
    PID pid_brk(600, 0, 50, 1500);
    PurePursuit pp;
    State s_pid = {0,0,0, 1.0, 0,0};

    ofstream f("comparison_results.csv");
    // Salviamo anche le coordinate di riferimento (Ref_X, Ref_Y) per il plot preciso
    f << "Time,Ref_Vx,Ref_X,Ref_Y,MPC_X,MPC_Y,MPC_Vx,MPC_Fx,MPC_Steer,PID_X,PID_Y,PID_Vx,PID_Fx,PID_Power\n";

    cout << "Simulazione Chicane START..." << endl;

    for(int i=0; i<steps; i++) {
        double t = i*dt;

        // --- MPC ---
        auto hor_mpc = get_horizon(s_mpc, track, op);
        double v_ref_mpc = hor_mpc.second; 
        Point current_ref = hor_mpc.first[0]; // Punto target attuale

        Eigen::VectorXd x0(6);
        x0 << s_mpc.x, s_mpc.y, s_mpc.theta, s_mpc.vx, s_mpc.vy, s_mpc.omega;
        mpc.updateDiscretization(x0, 0.0);
        
        pair<double,double> u_mpc;
        
        // Protezione contro errori solver
        try {
            u_mpc = mpc.compute(x0, hor_mpc.first, v_ref_mpc);
        } catch(...) {
            cout << "CRASH MPC al tempo " << t << endl;
            u_mpc = {0,0};
        }

        s_mpc = update_vehicle_physics(s_mpc, u_mpc.first, u_mpc.second, dt);

        // --- PID ---
        auto hor_pid = get_horizon(s_pid, track, op);
        double v_ref_pid = hor_pid.second;
        double steer_pid = pp.calc(s_pid, hor_pid.first);
        
        double fx_pid = 0;
        double err = v_ref_pid - s_pid.vx;
        if(err > 0.5) {
            fx_pid = pid_gas.calc(v_ref_pid, s_pid.vx, dt);
            pid_brk.reset();
        } else if(err < -0.5) {
            fx_pid = -pid_brk.calc(-v_ref_pid, -s_pid.vx, dt); 
            pid_gas.reset();
        }
        
        s_pid = update_vehicle_physics(s_pid, fx_pid, steer_pid, dt);

        // Log
        f << t << "," << v_ref_mpc << ","
          << current_ref.x << "," << current_ref.y << ","
          << s_mpc.x << "," << s_mpc.y << "," << s_mpc.vx << "," << u_mpc.first << "," << u_mpc.second << ","
          << s_pid.x << "," << s_pid.y << "," << s_pid.vx << "," << fx_pid << "," << (fx_pid*s_pid.vx) << "\n";

        if(s_mpc.x > 290) break;
    }
    f.close();
    cout << "Simulazione END. Dati in comparison_results.csv" << endl;
    return 0;
}