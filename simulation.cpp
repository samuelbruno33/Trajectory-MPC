#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "MPC.hpp"
#include "geometry.hpp"

using namespace std;

// ==========================================
// 1. STRUTTURE E UTILS
// ==========================================
struct State {
    double x, y, theta, vx, vy, omega;
};

double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// ==========================================
// 2. MODELLO VEICOLO (PLANT)
// ==========================================
State update_vehicle_physics(State s, double Fx, double delta, double dt) {
    double m = 320.0;
    double Iz = 98.03;
    double la = 0.792;
    double lb = 0.758;
    double B = 10.0, C = 1.9, D = 1.0; 

    // Saturazione Attuatori
    if (Fx > 1500) Fx = 1500;
    if (Fx < -1500) Fx = -1500;
    double max_steer = 24.0 * M_PI / 180.0;
    if (delta > max_steer) delta = max_steer;
    if (delta < -max_steer) delta = -max_steer;

    double vx_safe = std::max(s.vx, 1.0);
    double alpha_f = delta - atan2(s.vy + s.omega * la, vx_safe);
    double alpha_r = -atan2(s.vy - s.omega * lb, vx_safe);

    double Fyf = D * sin(C * atan(B * alpha_f)) * m * 9.81 * (lb / (la + lb)); 
    double Fyr = D * sin(C * atan(B * alpha_r)) * m * 9.81 * (la / (la + lb));

    double dot_vx = (Fx - Fyf * sin(delta)) / m + s.vy * s.omega;
    double dot_vy = (Fyf * cos(delta) + Fyr) / m - s.vx * s.omega;
    double dot_omega = (la * Fyf * cos(delta) - lb * Fyr) / Iz;

    State next_s = s;
    next_s.vx += dot_vx * dt;
    next_s.vy += dot_vy * dt;
    next_s.omega += dot_omega * dt;
    
    next_s.x += (s.vx * cos(s.theta) - s.vy * sin(s.theta)) * dt;
    next_s.y += (s.vx * sin(s.theta) + s.vy * cos(s.theta)) * dt;
    next_s.theta = normalize_angle(s.theta + s.omega * dt);

    return next_s;
}

// ==========================================
// 3. TRACK GENERATION (Rettilineo -> Sinusoide -> Rettilineo)
// ==========================================
struct TrackPoint {
    double x, y;
    double target_v; 
};

vector<TrackPoint> generate_track() {
    vector<TrackPoint> track;
    double ds = 0.5; 

    // A. Primo Rettilineo (0 -> 100m) - ACCELERAZIONE
    for (double x = 0; x < 100.0; x += ds) {
        track.push_back({x, 0.0, 25.0}); // Target 90 km/h
    }

    // B. Sinusoide (100 -> 200m) - GUIDATO LENTO
    // Ampiezza 4m, Frequenza adatta a FS
    for (double x = 100.0; x < 200.0; x += ds) {
        double y = 4.0 * sin((x - 100.0) / 15.0); 
        track.push_back({x, y, 12.0}); // Frenata a 43 km/h
    }

    // C. Secondo Rettilineo (200 -> 300m) - RI-ACCELERAZIONE
    for (double x = 200.0; x <= 300.0; x += ds) {
        track.push_back({x, 0.0, 25.0}); 
    }
    
    return track;
}

// Helper: Estrai orizzonte e riferimento velocità
pair<vector<Point>, double> get_local_horizon(State s, const vector<TrackPoint>& full_track, int horizon_len) {
    double min_dist = 1e9;
    int idx = 0;
    
    // Cerca il punto più vicino (Greedy search ottimizzata localmente si potrebbe fare, qui full scan)
    for (int i = 0; i < full_track.size(); i++) {
        double d = sqrt(pow(s.x - full_track[i].x, 2) + pow(s.y - full_track[i].y, 2));
        if (d < min_dist) { min_dist = d; idx = i; }
    }

    vector<Point> local_wp;
    // Riferimento velocità: prendiamo quello tra 1 secondo (circa 20m avanti) per anticipare la frenata
    // Se prendiamo quello attuale, frena troppo tardi!
    int lookahead_v_idx = std::min(idx + 15, (int)full_track.size() - 1);
    double target_v = full_track[lookahead_v_idx].target_v;

    for (int i = 0; i < horizon_len; i++) {
        int curr_idx = std::min(idx + i, (int)full_track.size() - 1);
        local_wp.push_back(Point(full_track[curr_idx].x, full_track[curr_idx].y));
    }

    return {local_wp, target_v};
}

// ==========================================
// 4. CONTROLLORI PID
// ==========================================
class PID {
    double kp, ki, kd, integral=0, prev_err=0, limit;
public:
    PID(double p, double i, double d, double lim) : kp(p), ki(i), kd(d), limit(lim) {}
    double update(double target, double current, double dt) {
        double err = target - current;
        integral += err * dt;
        double deriv = (err - prev_err) / dt;
        prev_err = err;
        
        // Anti-windup
        if(integral > 500) integral = 500; 
        if(integral < -500) integral = -500;

        double out = kp*err + ki*integral + kd*deriv;
        return std::max(-limit, std::min(out, limit));
    }
    void reset() { integral = 0; prev_err = 0; }
};

class PurePursuit {
public:
    double compute(State s, const vector<Point>& horizon) {
        double ld = 4.0; 
        for(auto& p : horizon) {
            double d = sqrt(pow(p.x - s.x, 2) + pow(p.y - s.y, 2));
            if(d > ld) {
                double alpha = atan2(p.y - s.y, p.x - s.x) - s.theta;
                return atan(2 * 1.55 * sin(alpha) / ld);
            }
        }
        return 0.0;
    }
};

// ==========================================
// 5. MAIN
// ==========================================
int main() {
    ifstream check("cornering_stiffness_vs_vertical_load.txt");
    if(!check.good()) cerr << "WARNING: File stiffness mancante!" << endl;

    auto track = generate_track();
    double dt = 0.05;
    double sim_time = 20.0; // Tempo sufficiente per fare 300m
    int steps = sim_time / dt;

    MPC mpc;
    // PID Tuning aggressivo per mostrare overshoot e violazioni
    PID pid_throttle(400, 2, 10, 1500); 
    PID pid_brake(600, 0, 20, 1500);
    PurePursuit pp;

    State s_mpc = {0,0,0, 1.0, 0,0}; 
    State s_pid = {0,0,0, 1.0, 0,0};

    ofstream file("comparison_results.csv");
    file << "Time,Ref_Vx,MPC_X,MPC_Y,MPC_Vx,MPC_Fx,MPC_Steer,MPC_LatErr,PID_X,PID_Y,PID_Vx,PID_Fx,PID_Power\n";

    cout << "Avvio simulazione Snake Track (Straight-Curve-Straight)..." << endl;

    for(int i=0; i<steps; i++) {
        double t = i*dt;

        // --- MPC ---
        auto horizon_mpc = get_local_horizon(s_mpc, track, 20);
        
        Eigen::VectorXd x0(6);
        x0 << s_mpc.x, s_mpc.y, s_mpc.theta, s_mpc.vx, s_mpc.vy, s_mpc.omega;
        mpc.updateDiscretization(x0, 0.0);
        
        pair<double, double> u_mpc = mpc.compute(x0, horizon_mpc.first, horizon_mpc.second);
        s_mpc = update_vehicle_physics(s_mpc, u_mpc.first, u_mpc.second, dt);

        // --- PID ---
        auto horizon_pid = get_local_horizon(s_pid, track, 20);
        double v_target_pid = horizon_pid.second;

        double steer_pp = pp.compute(s_pid, horizon_pid.first);
        
        double fx_pid = 0;
        double err_v = v_target_pid - s_pid.vx;
        
        if(err_v > 0.2) { // Deadzone piccola
            fx_pid = pid_throttle.update(v_target_pid, s_pid.vx, dt);
            pid_brake.reset();
        } else if (err_v < -0.2) {
            fx_pid = -pid_brake.update(-v_target_pid, -s_pid.vx, dt);
            pid_throttle.reset();
        }
        
        s_pid = update_vehicle_physics(s_pid, fx_pid, steer_pp, dt);

        // --- LOGGING ---
        // Calcolo LatErr approssimato (distanza Y dal centro pista nel tratto dritto o sinusoide)
        // Questo è solo per il plot veloce
        double ref_y_mpc = 0;
        if(s_mpc.x > 100 && s_mpc.x < 200) ref_y_mpc = 4.0 * sin((s_mpc.x - 100.0) / 15.0);
        double lat_err = s_mpc.y - ref_y_mpc;

        file << t << "," << horizon_mpc.second << ","
             << s_mpc.x << "," << s_mpc.y << "," << s_mpc.vx << "," << u_mpc.first << "," << u_mpc.second << ","
             << lat_err << "," 
             << s_pid.x << "," << s_pid.y << "," << s_pid.vx << "," << fx_pid << "," << (fx_pid * s_pid.vx) << "\n";
        
        // Stop se finita pista
        if(s_mpc.x > 300 && s_pid.x > 300) break;
    }

    file.close();
    cout << "Simulazione completata." << endl;
    return 0;
}