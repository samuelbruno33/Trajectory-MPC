#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <deque>
#include "MPC.hpp"
#include "geometry.hpp"

using namespace std;

// --- UTILS ---
double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// --- DELAY SIMULATOR (Per peggiorare il PID) ---
class ActuatorDelay {
    std::deque<double> buffer;
    int delay_steps;
public:
    ActuatorDelay(int steps) : delay_steps(steps) {
        for(int i=0; i<steps; i++) buffer.push_back(0.0);
    }
    double get_delayed_value(double new_input) {
        if (delay_steps <= 0) return new_input;
        buffer.push_back(new_input);
        double val = buffer.front();
        buffer.pop_front();
        return val;
    }
};

// --- PLANT ---
struct State { double x, y, theta, vx, vy, omega; };

State update_vehicle_physics(State s, double Fx, double delta, double dt) {
    double m = 320.0; double Iz = 98.03; double la = 0.792; double lb = 0.758;
    double B = 10.0, C = 1.9, D = 1.0; 

    // Saturazione
    if (Fx > 1500) Fx = 1500; if (Fx < -1500) Fx = -1500;
    double max_steer = 0.42; 
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

// --- TRACK (Sinusoide) ---
struct TrackPoint { double x, y, v_target; };

vector<TrackPoint> generate_track() {
    vector<TrackPoint> t;
    double ds = 0.5; 

    // Rettilineo
    for (double x = 0; x < 40.0; x += ds) {
        t.push_back({x, 0.0, 25.0}); 
    }
    // Sinusoide (Curvoni)
    double x_start_sine = 40.0;
    double length_sine = 200.0;
    
    for (double x = 40.0; x < 240.0; x += ds) {
        double angle = ((x - x_start_sine) / length_sine) * 2.0 * M_PI;
        double y = 6.0 * sin(angle); 
        t.push_back({x, y, 15.0}); 
    }
    // Rettilineo Finale
    for (double x = 240.0; x <= 300.0; x += ds) {
        t.push_back({x, 0.0, 25.0}); 
    }
    return t;
}

// --- HORIZON (Time-Based) ---
pair<vector<Point>, double> get_horizon(State s, const vector<TrackPoint>& track, int op, double dt) {
    double min_d = 1e9; int current_idx = 0;
    for(int i=0; i<track.size(); i++) {
        double d = sqrt(pow(s.x - track[i].x, 2) + pow(s.y - track[i].y, 2));
        if(d < min_d) { min_d = d; current_idx = i; }
    }
    
    vector<Point> wp;
    double current_v_ref = track[current_idx].v_target;
    double track_res = 0.5; 
    double proj_speed = std::max(s.vx, 5.0);

    for(int i=0; i<op; i++) {
        double future_dist = proj_speed * dt * (i + 1);
        int idx_offset = static_cast<int>(future_dist / track_res);
        int target_idx = std::min(current_idx + idx_offset, (int)track.size()-1);
        wp.push_back(Point(track[target_idx].x, track[target_idx].y));
    }
    return {wp, current_v_ref};
}

// --- PID (Setup ORIGINALE GITHUB) ---
class PID {
    double kp, ki, kd, I=0, prev=0, limit;
public:
    PID(double p, double i, double d, double lim) : kp(p), ki(i), kd(d), limit(lim) {}
    double calc(double ref, double meas, double dt) {
        double err = ref - meas;
        I += err*dt;
        if(I > 300) I=300; if(I<-300) I=-300; 
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
        double ld = 5.0; // Valore GitHub (Standard)
        
        for(auto& p : path) {
            double dist = sqrt(pow(p.x-s.x, 2) + pow(p.y-s.y, 2));
            if(dist > ld) {
                double alpha = atan2(p.y - s.y, p.x - s.x) - s.theta;
                alpha = normalize_angle(alpha);
                return atan(2 * 0.792 * sin(alpha) / ld);
            }
        }
        return 0;
    }
};

int main() {
    ifstream check("cornering_stiffness_vs_vertical_load.txt");
    if(!check.good()) { cerr << "ERRORE: File Dati Mancante!" << endl; return -1; }

    auto track = generate_track();
    int op = 40; 
    double dt = 0.05;
    int steps = 600; 

    MPC mpc; 
    State s_mpc = {0,0,0, 1.0, 0,0};
    
    // --- PID CONFIG (RIPRESA DA GITHUB) ---
    // Gas: P=400, I=5, D=20 (Molto più morbido del mio 1000)
    // Brake: P=600, I=0, D=50
    PID pid_gas(400, 5, 20, 1500); 
    PID pid_brk(600, 0, 50, 1500);
    PurePursuit pp;
    State s_pid = {0,0,0, 1.0, 0,0};

    // --- RITARDI ---
    // Ritardo sterzo a 2 steps (sufficiente per creare oscillazione e ritardo in curva)
    ActuatorDelay pid_steer_delay(2); 
    ActuatorDelay pid_gas_delay(0);

    ofstream f("comparison_results.csv");
    f << "Time,Ref_Vx,Ref_X,Ref_Y,MPC_X,MPC_Y,MPC_Vx,MPC_Fx,MPC_Steer,PID_X,PID_Y,PID_Vx,PID_Fx,PID_Power\n";

    cout << "Simulazione START..." << endl;

    for(int i=0; i<steps; i++) {
        double t = i*dt;

        // --- MPC ---
        auto hor_mpc = get_horizon(s_mpc, track, op, dt); 
        double v_ref_mpc = hor_mpc.second; 
        
        Eigen::VectorXd x0(6);
        x0 << s_mpc.x, s_mpc.y, s_mpc.theta, s_mpc.vx, s_mpc.vy, s_mpc.omega;
        mpc.updateDiscretization(x0, 0.0);
        
        pair<double,double> u_mpc;
        try { u_mpc = mpc.compute(x0, hor_mpc.first, v_ref_mpc); } catch(...) { u_mpc = {0,0}; }
        s_mpc = update_vehicle_physics(s_mpc, u_mpc.first, u_mpc.second, dt);

        // --- PID ---
        auto hor_pid = get_horizon(s_pid, track, op, dt);
        double steer_pid_req = pp.calc(s_pid, hor_pid.first);
        double fx_pid_req = 0;
        double err_v = hor_pid.second - s_pid.vx;
        
        if(err_v > 0.2) {
            fx_pid_req = pid_gas.calc(hor_pid.second, s_pid.vx, dt);
            pid_brk.reset();
        } else if(err_v < -0.2) {
            fx_pid_req = -pid_brk.calc(-hor_pid.second, -s_pid.vx, dt); 
            pid_gas.reset();
        }

        double real_steer_pid = pid_steer_delay.get_delayed_value(steer_pid_req);
        double real_fx_pid = pid_gas_delay.get_delayed_value(fx_pid_req);
        s_pid = update_vehicle_physics(s_pid, real_fx_pid, real_steer_pid, dt);

        f << t << "," << v_ref_mpc << ","
          << hor_mpc.first[0].x << "," << hor_mpc.first[0].y << ","
          << s_mpc.x << "," << s_mpc.y << "," << s_mpc.vx << "," << u_mpc.first << "," << u_mpc.second << ","
          << s_pid.x << "," << s_pid.y << "," << s_pid.vx << "," << real_fx_pid << "," << (real_fx_pid*s_pid.vx) << "\n";

        if(s_mpc.x > 295) break;
    }
    f.close();
    cout << "Simulazione END." << endl;
    return 0;
}