#include "MPC.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// Header fondamentale per usare .exp() sulle matrici Eigen
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;

MPC::MPC() : solver_initialized(false) {
    // --- 1. INIZIALIZZAZIONE FORZATA VARIABILI ---
    nx = 6;  // [X, Y, Theta, Vx, Vy, Omega]
    nu = 2;  // [Fx, Delta]
    op = 40; // Orizzonte Predizione
    oc = 10; // Orizzonte Controllo
    dt = 0.05; 

    // Parametri Veicolo
    la = 0.792; lb = 0.758; Iz = 98.03; m = 320; Ycm = 0.362;
    P_max = 30000.0; F_peak = 1500.0;

    // --- 2. SETUP MATRICI PESI (TUNING OTTIMIZZATO) ---
    Q = Eigen::MatrixXd::Zero(nx, nx);

    // ✅ TUNING MIGLIORATO (TESI-READY):
    // X: 80 → Aumentato per ridurre "taglio curva" finale (250-300m)
    // Y: 100 → Precisione laterale mantenuta
    // Theta: 60 → Aumentato per allineamento in curva
    // Vx: 250 → Tracking velocità prioritario
    // Vy, Omega: Bassi per permettere derive controllate
    Q.diagonal() << 80, 100, 60, 250, 10, 1; 

    R = Eigen::MatrixXd::Identity(nu, nu);
    // Costo input ridotto per risposta aggressiva
    R.diagonal() << 5e-5, 0.05; 

    R_delta = Eigen::MatrixXd::Identity(nu, nu); 
    // Penalità variazione: permette scatti rapidi ma evita chattering
    R_delta.diagonal() << 1e-5, 0.8; 

    // --- 3. SETUP LIMITI ---
    umin = Eigen::VectorXd::Zero(nu);
    umax = Eigen::VectorXd::Zero(nu);
    
    umin[IDX_DELTA] = -24.0 * M_PI / 180.0; 
    umax[IDX_DELTA] =  24.0 * M_PI / 180.0; 
    umin[IDX_FX]    = -F_peak; 
    umax[IDX_FX]    =  F_peak; 

    xmin = Eigen::VectorXd::Constant(nx, -OsqpEigen::INFTY);
    xmax = Eigen::VectorXd::Constant(nx, OsqpEigen::INFTY);
    u_prev = Eigen::VectorXd::Zero(oc * nu);
    
    numeri = load_data();
}

void MPC::updateDiscretization(const Eigen::VectorXd& x, double acc_y_measured) {
    A = Eigen::MatrixXd::Zero(nx, nx);
    B = Eigen::MatrixXd::Zero(nx, nu);

    double theta = x(IDX_THETA);
    double vx = std::max(x(IDX_VX), 1.0); 
    double vy = x(IDX_VY);
    double omega = x(IDX_OMEGA);
    double sin_th = sin(theta);
    double cos_th = cos(theta);

    double acc_x_est = 0.0; 
    std::pair<double, double> K_values = load_transfer(acc_x_est, numeri);
    double Ka = K_values.first;  
    double Kp = K_values.second; 

    // --- MATRICE A ---
    A(IDX_X, IDX_THETA) = -vx * sin_th - vy * cos_th;
    A(IDX_X, IDX_VX)    = cos_th;
    A(IDX_X, IDX_VY)    = -sin_th;

    A(IDX_Y, IDX_THETA) = vx * cos_th - vy * sin_th;
    A(IDX_Y, IDX_VX)    = sin_th;
    A(IDX_Y, IDX_VY)    = cos_th;

    A(IDX_THETA, IDX_OMEGA) = 1.0;

    A(IDX_VX, IDX_VY)    = -omega; 
    A(IDX_VX, IDX_OMEGA) = -vy;

    double den = m * vx;
    double C_vy = -(Ka + Kp) / den; 
    double C_omega = (Kp * lb - Ka * la) / den - vx; 
    
    A(IDX_VY, IDX_VY)    = C_vy;
    A(IDX_VY, IDX_OMEGA) = C_omega;
    A(IDX_VY, IDX_VX)    = -omega; 

    double den_rot = Iz * vx;
    double D_vy = (Kp * lb - Ka * la) / den_rot;
    double D_omega = -(Ka * la * la + Kp * lb * lb) / den_rot;

    A(IDX_OMEGA, IDX_VY)    = D_vy;
    A(IDX_OMEGA, IDX_OMEGA) = D_omega;

    // --- MATRICE B ---
    B(IDX_VX, IDX_FX) = 1.0 / m;
    B(IDX_VY, IDX_DELTA)    = Ka / m;
    B(IDX_OMEGA, IDX_DELTA) = (Ka * la) / Iz;

    // --- DISCRETIZZAZIONE ZOH ---
    Eigen::MatrixXd M(nx + nu, nx + nu);
    M.setZero();
    M.topLeftCorner(nx, nx) = A * dt;
    M.topRightCorner(nx, nu) = B * dt;
    
    Eigen::MatrixXd expM = M.exp();
    Ad = expM.topLeftCorner(nx, nx);
    Bd = expM.topRightCorner(nx, nu);
}

std::pair<double, double> MPC::compute(const Eigen::VectorXd& x0, const vector<Point>& waypoints, double v_ref) {
    // 1. Vincolo Potenza Dinamico
    double current_vx = std::max(x0(IDX_VX), 1.0);
    double max_force_power = P_max / current_vx;
    double current_fx_max = std::min(F_peak, max_force_power);

    // 2. Matrici Blocchi
    Eigen::MatrixXd Q_blk = Eigen::MatrixXd::Zero(op * nx, op * nx);
    Eigen::MatrixXd R_blk = Eigen::MatrixXd::Zero(oc * nu, oc * nu);
    Eigen::MatrixXd Rd_blk = Eigen::MatrixXd::Zero((oc - 1) * nu, (oc - 1) * nu);

    for (int i = 0; i < op; ++i) Q_blk.block(i * nx, i * nx, nx, nx) = Q; 
    for (int i = 0; i < oc; ++i) R_blk.block(i * nu, i * nu, nu, nu) = R;
    for (int i = 0; i < oc - 1; ++i) Rd_blk.block(i * nu, i * nu, nu, nu) = R_delta;

    // 3. Matrici Predittive
    Eigen::MatrixXd Sx = Eigen::MatrixXd::Zero(op * nx, nx); 
    Eigen::MatrixXd Su = Eigen::MatrixXd::Zero(op * nx, oc * nu);
    Eigen::MatrixXd A_power = Eigen::MatrixXd::Identity(nx, nx);
    
    for (int i = 0; i < op; ++i) {
        A_power *= Ad;
        Sx.block(i * nx, 0, nx, nx) = A_power;
        if (i < oc) {
             for (int j = 0; j <= i; ++j) {
                 Eigen::MatrixXd term = Eigen::MatrixXd::Identity(nx, nx);
                 for(int k=0; k < (i-j); ++k) term *= Ad;
                 Su.block(i*nx, j*nu, nx, nu) = term * Bd;
             }
        }
    }

    // 4. ✅ COSTRUZIONE RIFERIMENTO MIGLIORATA
    Eigen::VectorXd x_ref_big = Eigen::VectorXd::Zero(op * nx);
    
    // Calcolo theta_ref tramite media mobile (riduce rumore da waypoints ravvicinati)
    for (int i = 0; i < op; ++i) {
        int idx = std::min(i, (int)waypoints.size() - 1);
        
        x_ref_big(i * nx + IDX_X) = waypoints[idx].x;
        x_ref_big(i * nx + IDX_Y) = waypoints[idx].y;
        
        // ✅ Theta smoothing: usa differenza tra punto i+5 e i-5 (se disponibili)
        int idx_ahead = std::min(idx + 5, (int)waypoints.size() - 1);
        int idx_behind = std::max(idx - 5, 0);
        
        double dx = waypoints[idx_ahead].x - waypoints[idx_behind].x;
        double dy = waypoints[idx_ahead].y - waypoints[idx_behind].y;
        double ref_theta = 0.0;
        
        if(sqrt(dx*dx + dy*dy) > 0.5) {
            ref_theta = atan2(dy, dx);
        } else if (i > 0) {
            ref_theta = x_ref_big((i-1) * nx + IDX_THETA);
        } else {
            ref_theta = x0(IDX_THETA);
        }
        
        x_ref_big(i * nx + IDX_THETA) = ref_theta;
        x_ref_big(i * nx + IDX_VX) = v_ref; 
    }

    // 5. QP Setup
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero((oc - 1) * nu, oc * nu);
    for (int i = 0; i < oc - 1; ++i) {
        D.block(i * nu, i * nu, nu, nu) = -Eigen::MatrixXd::Identity(nu, nu); 
        D.block(i * nu, (i + 1) * nu, nu, nu) = Eigen::MatrixXd::Identity(nu, nu); 
    }

    Eigen::MatrixXd H = Su.transpose() * Q_blk * Su + R_blk + D.transpose() * Rd_blk * D;
    Eigen::VectorXd x_pred_free = Sx * x0; 
    Eigen::VectorXd error_term = (x_pred_free - x_ref_big).transpose() * Q_blk * Su;
    Eigen::VectorXd slew_term = Eigen::VectorXd::Zero(oc * nu);
    slew_term.head(nu) = -u_prev.head(nu).transpose() * R_delta; 
    Eigen::VectorXd q = error_term + slew_term;

    // 6. Vincoli (Safety Checks Inclusi)
    int n_constraints = oc * nu + (oc - 1) * nu; 
    Eigen::MatrixXd A_total = Eigen::MatrixXd::Zero(n_constraints, oc * nu);
    Eigen::VectorXd lb_total = Eigen::VectorXd::Zero(n_constraints);
    Eigen::VectorXd ub_total = Eigen::VectorXd::Zero(n_constraints);

    // Vincoli Input
    A_total.topRows(oc * nu) = Eigen::MatrixXd::Identity(oc * nu, oc * nu);
    for(int i=0; i<oc; i++) {
        lb_total(i*nu + IDX_FX) = -current_fx_max; 
        ub_total(i*nu + IDX_FX) = current_fx_max;
        lb_total(i*nu + IDX_DELTA) = umin[IDX_DELTA];
        ub_total(i*nu + IDX_DELTA) = umax[IDX_DELTA];
        
        // Safety swap
        if (lb_total(i*nu + IDX_FX) > ub_total(i*nu + IDX_FX)) {
             double t = lb_total(i*nu + IDX_FX);
             lb_total(i*nu + IDX_FX) = ub_total(i*nu + IDX_FX);
             ub_total(i*nu + IDX_FX) = t;
        }
    }

    // Vincoli Slew Rate (Realistici per attuatori reali)
    double max_dFx = 8000.0 * dt;      // ✅ 8kN/s (realistico per motori EV)
    double max_dDelta = 50.0 * dt;     // ✅ 50 rad/s (realistico per sterzo EPS)
    
    A_total.block(oc*nu, 0, (oc-1)*nu, oc*nu) = D;
    for(int i=0; i < oc-1; i++) {
        lb_total(oc*nu + i*nu + IDX_FX) = -max_dFx;
        ub_total(oc*nu + i*nu + IDX_FX) = max_dFx;
        lb_total(oc*nu + i*nu + IDX_DELTA) = -max_dDelta;
        ub_total(oc*nu + i*nu + IDX_DELTA) = max_dDelta;
    }

    Eigen::SparseMatrix<double> H_s = H.sparseView();
    Eigen::SparseMatrix<double> A_s = A_total.sparseView();

    // 7. Solver Lifecycle
    if (!solver_initialized) {
        // CLEANUP CRUCIALE
        solver.data()->clearHessianMatrix();
        solver.data()->clearLinearConstraintsMatrix();

        solver.settings()->setWarmStart(true);
        solver.settings()->setVerbosity(false);
        solver.settings()->setScaling(0); 
        
        solver.data()->setNumberOfVariables(oc * nu);
        solver.data()->setNumberOfConstraints(n_constraints);
        
        if (!solver.data()->setHessianMatrix(H_s)) return {0.0, 0.0};
        if (!solver.data()->setGradient(q)) return {0.0, 0.0};
        if (!solver.data()->setLinearConstraintsMatrix(A_s)) return {0.0, 0.0};
        if (!solver.data()->setLowerBound(lb_total)) return {0.0, 0.0};
        if (!solver.data()->setUpperBound(ub_total)) return {0.0, 0.0};
        
        if (!solver.initSolver()) return {0.0, 0.0};
        solver_initialized = true;
    } else {
        solver.updateHessianMatrix(H_s);
        solver.updateGradient(q);
        solver.updateBounds(lb_total, ub_total);
    }

    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
        solver_initialized = false;
        return {0.0, 0.0}; 
    }

    Eigen::VectorXd u_opt = solver.getSolution();
    u_prev = u_opt; 
    return { u_opt(0), u_opt(1) };
}

vector<vector<double>> MPC::load_data() {
    ifstream file("cornering_stiffness_vs_vertical_load.txt");
    if (!file) return {{0.0, 0.0}, {2000.0, 100000.0}}; 
    vector<vector<double>> nums;
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        double val1, val2;
        if (ss >> val1 >> val2) nums.push_back({val1, val2});
    }
    file.close();
    return nums;
}

double interpolate(double Fx, const vector<vector<double>> &table) {
    if (table.empty()) return 0.0;
    if (Fx <= table.front()[0]) return table.front()[1];
    if (Fx >= table.back()[0]) return table.back()[1];
    for (size_t i = 1; i < table.size(); ++i) {
        if (Fx >= table[i - 1][0] && Fx <= table[i][0]) {
            double F1 = table[i - 1][0];
            double F2 = table[i][0];
            double C1 = table[i - 1][1];
            double C2 = table[i][1];
            return C1 + (((Fx - F1) / (F2 - F1)) * (C2 - C1));
        }
    }
    return 0.0;
}

pair<double,double> MPC::load_transfer(double acceleration_x, const vector<vector<double>> &numeri){
    double L = la + lb;
    double weight_term = m * 9.81;
    double transfer_term = m * acceleration_x * Ycm;
    double Fant = (weight_term * lb - transfer_term) / L;
    double Frear = (weight_term * la + transfer_term) / L;
    double Ka = interpolate(Fant, numeri);
    double Kp = interpolate(Frear, numeri);
    return make_pair(Ka, Kp);
}
