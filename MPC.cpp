#include "MPC.hpp"
#include <algorithm>
#include <cmath>

// Header fondamentale per usare .exp() sulle matrici Eigen
#include <unsupported/Eigen/MatrixFunctions>

MPC::MPC() : solver_initialized(false) {
    // --- Inizializzazione Pesi (Tuning) ---
    
    // Q: Matrice dei pesi sugli stati [X, Y, Theta, Vx, Vy, Omega]
    Q = Eigen::MatrixXd::Zero(nx, nx);
    // Tuning preliminare:
    // - Vx: peso alto (50) per garantire il tracking della velocità
    // - Y: peso medio (20) per minimizzare l'errore laterale
    // - Vy: peso per la stabilità (10)
    Q.diagonal() << 20, 20, 10, 50, 10, 1; 

    // R: Matrice dei pesi sugli input [Fx, Delta]
    R = Eigen::MatrixXd::Identity(nu, nu);
    // - Fx: peso basso (1e-4) per non limitare troppo l'accelerazione
    // - Delta: peso alto (1.0) per penalizzare l'uso eccessivo dello sterzo
    R.diagonal() << 1e-4, 1.0; 

    // R_delta: Matrice dei pesi sulla variazione degli input (Slew Rate)
    R_delta = Eigen::MatrixXd::Identity(nu, nu); 
    // - Delta Fx: peso basso
    // - Delta Sterzo: peso molto alto (100) per evitare oscillazioni/chattering
    R_delta.diagonal() << 1e-3, 100.0; 

    // --- Definizione Limiti Hard (Costanti) ---
    umin = Eigen::VectorXd::Zero(nu);
    umax = Eigen::VectorXd::Zero(nu);
    
    // Vincoli Sterzo (fisici)
    umin[IDX_DELTA] = -24.0 * M_PI / 180.0; // -24 gradi in radianti
    umax[IDX_DELTA] = 24.0 * M_PI / 180.0; // +24 gradi in radianti
    
    // Vincoli Forza (inizializzazione al picco, poi limitati dinamicamente in compute)
    umin[IDX_FX] = -F_peak; // Massima frenata
    umax[IDX_FX] = F_peak; // Massima trazione

    // Inizializzazione vettori stato limiti (per ora infiniti, gestiti nei vincoli lineari)
    xmin = Eigen::VectorXd::Constant(nx, -OsqpEigen::INFTY);
    xmax = Eigen::VectorXd::Constant(nx, OsqpEigen::INFTY);

    // Inizializzazione vettore warm start
    u_prev = Eigen::VectorXd::Zero(oc * nu);
    
    // Caricamento dati pneumatici
    numeri = load_data();
}

void MPC::updateDiscretization(const Eigen::VectorXd& x, double acc_y_measured) {
    A = Eigen::MatrixXd::Zero(nx, nx);
    B = Eigen::MatrixXd::Zero(nx, nu);

    // Estrazione stato corrente
    double theta = x(IDX_THETA);
    double vx = std::max(x(IDX_VX), 1.0); // Clamp a 1 m/s per evitare divisioni per zero
    double vy = x(IDX_VY);
    double omega = x(IDX_OMEGA);

    double sin_th = sin(theta);
    double cos_th = cos(theta);

    // Calcolo Cornering Stiffness basata sul trasferimento di carico
    // Nota: usiamo acc_x nullo per semplificazione o stimato se disponibile
    double acc_x_est = 0.0; 
    std::pair<double, double> K_values = load_transfer(acc_x_est, numeri);
    double Ka = K_values.first;  // Stiffness anteriore (K1)
    double Kp = K_values.second; // Stiffness posteriore (K2)

    // --- COSTRUZIONE MATRICE A (Jacobiana 6x6) ---
    // Linearizzazione del modello a bicicletta dinamico accoppiato

    // 1. Equazioni Posizione (Cinematica Rotata)
    // dX = vx*cos(th) - vy*sin(th)
    A(IDX_X, IDX_THETA) = -vx * sin_th - vy * cos_th;
    A(IDX_X, IDX_VX)    = cos_th;
    A(IDX_X, IDX_VY)    = -sin_th;

    // dY = vx*sin(th) + vy*cos(th)
    A(IDX_Y, IDX_THETA) = vx * cos_th - vy * sin_th;
    A(IDX_Y, IDX_VX)    = sin_th;
    A(IDX_Y, IDX_VY)    = cos_th;

    // 2. Equazione Imbardata
    A(IDX_THETA, IDX_OMEGA) = 1.0;

    // 3. Dinamica Longitudinale (Accoppiata)
    // dVx = Fx/m - vy*omega
    A(IDX_VX, IDX_VY)    = -omega; // Termine Coriolis
    A(IDX_VX, IDX_OMEGA) = -vy;

    // 4. Dinamica Laterale
    // dVy = (Fyf + Fyr)/m - vx*omega
    double den = m * vx;
    double C_vy = -(Ka + Kp) / den; // Coeff per vy
    // Coeff per omega: include termine dinamico e termine cinematico (-vx)
    double C_omega = (Kp * lb - Ka * la) / den - vx; 
    
    A(IDX_VY, IDX_VY)    = C_vy;
    A(IDX_VY, IDX_OMEGA) = C_omega;
    // Approssimazione: trascuriamo la derivata parziale rispetto a Vx (1/v^2) per stabilità numerica
    A(IDX_VY, IDX_VX)    = -omega; 

    // 5. Dinamica Rotazionale
    // dOmega = (a*Fyf - b*Fyr)/Iz
    double den_rot = Iz * vx;
    double D_vy = (Kp * lb - Ka * la) / den_rot;
    double D_omega = -(Ka * la * la + Kp * lb * lb) / den_rot;

    A(IDX_OMEGA, IDX_VY)    = D_vy;
    A(IDX_OMEGA, IDX_OMEGA) = D_omega;

    // --- COSTRUZIONE MATRICE B (Jacobiana 6x2) ---
    
    // Ingresso Forza Longitudinale (Fx)
    B(IDX_VX, IDX_FX) = 1.0 / m;

    // Ingresso Sterzo (Delta) -> influenza Vy e Omega tramite Fyf
    B(IDX_VY, IDX_DELTA)    = Ka / m;
    B(IDX_OMEGA, IDX_DELTA) = (Ka * la) / Iz;

    // --- DISCRETIZZAZIONE ZOH (Esatta tramite esponenziale di matrice) ---
    Eigen::MatrixXd M(nx + nu, nx + nu);
    M.setZero();
    M.topLeftCorner(nx, nx) = A * dt;
    M.topRightCorner(nx, nu) = B * dt;
    
    Eigen::MatrixXd expM = M.exp();

    // Estrazione matrici discrete
    Ad = expM.topLeftCorner(nx, nx);
    Bd = expM.topRightCorner(nx, nu);
}

std::pair<double, double> MPC::compute(const Eigen::VectorXd& x0, const vector<Point>& waypoints, double v_ref) {
        
    // 1. Calcolo Vincolo Potenza Dinamico
    // Fx_max <= P_max / v_current
    double current_vx = std::max(x0(IDX_VX), 1.0);
    double max_force_power = P_max / current_vx;
    double current_fx_max = std::min(F_peak, max_force_power);

    // 2. Costruzione Matrici a Blocchi (Prediction Matrices)
    Eigen::MatrixXd Q_blk = Eigen::MatrixXd::Zero(op * nx, op * nx);
    Eigen::MatrixXd R_blk = Eigen::MatrixXd::Zero(oc * nu, oc * nu);
    Eigen::MatrixXd Rd_blk = Eigen::MatrixXd::Zero((oc - 1) * nu, (oc - 1) * nu);

    for (int i = 0; i < op; ++i) 
        Q_blk.block(i * nx, i * nx, nx, nx) = Q; 
    for (int i = 0; i < oc; ++i) 
        R_blk.block(i * nu, i * nu, nu, nu) = R;
    for (int i = 0; i < oc - 1; ++i) 
        Rd_blk.block(i * nu, i * nu, nu, nu) = R_delta;

    // 3. Matrici Predittive Sx (Evoluzione Libera) e Su (Convoluzione)
    Eigen::MatrixXd Sx = Eigen::MatrixXd::Zero(op * nx, nx); 
    Eigen::MatrixXd Su = Eigen::MatrixXd::Zero(op * nx, oc * nu);

    Eigen::MatrixXd A_power = Eigen::MatrixXd::Identity(nx, nx);
    for (int i = 0; i < op; ++i) {
        A_power *= Ad;
        Sx.block(i * nx, 0, nx, nx) = A_power; // A^1, A^2... A^op
        
        // Riempimento Su (triangolare inferiore a blocchi)
        if (i < oc) {
             // Calcolo somma di convoluzione
             for (int j = 0; j <= i; ++j) {
                 // Termine A^(i-j)*B
                 Eigen::MatrixXd term = Eigen::MatrixXd::Identity(nx, nx);
                 for(int k=0; k < (i-j); ++k) term *= Ad; // Potenza inefficiente ma chiara
                 Su.block(i*nx, j*nu, nx, nu) = term * Bd;
             }
        } else {
             // Se op > oc, continuiamo a propagare l'ultimo input
             // (Qui semplificato assumendo Su dimensionato correttamente per oc)
             // Per l'implementazione base, assumiamo che l'influenza continui
             // ma per ora manteniamo la struttura standard.
        }
    }

    // 4. Costruzione Vettore Riferimento (Trajectory & Speed Tracking)
    Eigen::VectorXd x_ref_big = Eigen::VectorXd::Zero(op * nx);
    for (int i = 0; i < op; ++i) {
        int idx = std::min(i, (int)waypoints.size() - 1);
        x_ref_big(i * nx + IDX_X) = waypoints[idx].x;
        x_ref_big(i * nx + IDX_Y) = waypoints[idx].y;
        // theta ref approssimato o calcolato dalla geometria waypoints (qui 0 o interpolato)
        x_ref_big(i * nx + IDX_VX) = v_ref; // Tracking Velocità Fondamentale
        // Vy e Omega rif = 0 (per stabilità)
    }

    // 5. Costruzione Problema QP: 1/2 U'HU + q'U
    // H = Su'QSu + R + D'RdD
    
    // Matrice D per Slew Rate (Delta U)
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero((oc - 1) * nu, oc * nu);
    for (int i = 0; i < oc - 1; ++i) {
        D.block(i * nu, i * nu, nu, nu) = -Eigen::MatrixXd::Identity(nu, nu); 
        D.block(i * nu, (i + 1) * nu, nu, nu) = Eigen::MatrixXd::Identity(nu, nu); 
    }

    Eigen::MatrixXd H = Su.transpose() * Q_blk * Su + R_blk + D.transpose() * Rd_blk * D;
    
    // Termine lineare q
    Eigen::VectorXd x_pred_free = Sx * x0; // Evoluzione libera
    Eigen::VectorXd error_term = (x_pred_free - x_ref_big).transpose() * Q_blk * Su;
    
    // Termine relativo a u_prev nel costo slew rate
    // Espansione: (u0 - u_prev)^T R_delta (u0 - u_prev) -> -2 u_prev^T R_delta u0
    // Aggiungiamo il contributo al gradiente per il primo elemento u0
    Eigen::VectorXd slew_term = Eigen::VectorXd::Zero(oc * nu);
    slew_term.head(nu) = -u_prev.head(nu).transpose() * R_delta; // Semplificazione primo passo
    
    // Nota: La formula corretta completa per q col termine D matrix è complessa, 
    // spesso si approssima o si usa il termine separato. Qui usiamo la versione base.
    Eigen::VectorXd q = error_term + slew_term; // + termine dovuto a u_prev

    // 6. Gestione Vincoli (A_c U <= bounds)
    // Vincoli considerati: 
    // a. Limiti Input (Hard)
    // b. Slew Rate (Delta U)
    // c. Corridoio Posizione (Y tracking) - Qui semplificato come vincoli su stati
    
    // Matrice Vincoli Totale A_total
    // Struttura:
    // [ I_u  ]  (Input constraints)
    // [ D    ]  (Slew Rate constraints)
    
    int n_constraints = oc * nu + (oc - 1) * nu; 
    Eigen::MatrixXd A_total = Eigen::MatrixXd::Zero(n_constraints, oc * nu);
    Eigen::VectorXd lb_total = Eigen::VectorXd::Zero(n_constraints);
    Eigen::VectorXd ub_total = Eigen::VectorXd::Zero(n_constraints);

    // a. Vincoli Input (Input Constraints)
    A_total.topRows(oc * nu) = Eigen::MatrixXd::Identity(oc * nu, oc * nu);
    
    for(int i=0; i<oc; i++) {
        // Fx bounds (Dinamici)
        lb_total(i*nu + IDX_FX) = -current_fx_max; 
        ub_total(i*nu + IDX_FX) = current_fx_max;
        
        // Delta bounds (Statici)
        lb_total(i*nu + IDX_DELTA) = umin[IDX_DELTA];
        ub_total(i*nu + IDX_DELTA) = umax[IDX_DELTA];
    }

    // b. Vincoli Slew Rate (Delta U constraints)
    // Max variazione per step (dt = 0.05s)
    double max_dFx = 3000.0 * dt; // N/s * s = N max variation
    double max_dDelta = 5.0 * dt; // rad/s * s
    
    A_total.block(oc*nu, 0, (oc-1)*nu, oc*nu) = D;
    
    for(int i=0; i < oc-1; i++) {
        lb_total(oc*nu + i*nu + IDX_FX) = -max_dFx;
        ub_total(oc*nu + i*nu + IDX_FX) = max_dFx;
        
        lb_total(oc*nu + i*nu + IDX_DELTA) = -max_dDelta;
        ub_total(oc*nu + i*nu + IDX_DELTA) = max_dDelta;
    }

    // Conversione Sparse per OSQP
    Eigen::SparseMatrix<double> H_s = H.sparseView();
    Eigen::SparseMatrix<double> A_s = A_total.sparseView();

    // 7. Solver Setup & Solve
    if (!solver_initialized) {
        solver.settings()->setWarmStart(true);
        solver.settings()->setVerbosity(true); // Silenzioso per real-time, TODO: per test tenere a true, cambiare in RT
        
        solver.data()->setNumberOfVariables(oc * nu);
        solver.data()->setNumberOfConstraints(n_constraints);
        
        if (!solver.data()->setHessianMatrix(H_s)) return {0.0, 0.0};
        if (!solver.data()->setGradient(q)) return {0.0, 0.0};
        if (!solver.data()->setLinearConstraintsMatrix(A_s)) return {0.0, 0.0};
        if (!solver.data()->setLowerBound(lb_total)) return {0.0, 0.0};
        if (!solver.data()->setUpperBound(ub_total)) return {0.0, 0.0};

        if (!solver.initSolver()) {
            std::cerr << "OSQP solver init failed\n";
            return {0.0, 0.0};
        }
        solver_initialized = true;
    } else {
        if (!solver.updateHessianMatrix(H_s)) return {0.0, 0.0};
        if (!solver.updateGradient(q)) return {0.0, 0.0};
        if (!solver.updateBounds(lb_total, ub_total)) return {0.0, 0.0};
    }

    // Risoluzione
    if (!solver.solve()) { 
        std::cerr << "OSQP: Risoluzione fallita\n";
        return {0.0, 0.0}; 
    }

    // Estrazione Soluzione
    Eigen::VectorXd u_opt = solver.getSolution();
    
    // Aggiornamento Warm Start
    u_prev = u_opt; 

    // Return primo input ottimale: Fx, Delta
    return { u_opt(0), u_opt(1) };
}

vector<vector<double>> MPC::load_data() {
    ifstream file("cornering_stiffness_vs_vertical_load.txt");
    if (!file) {
        cerr << "Impossibile aprire il file cornering stiffness!" << endl;
        // Valori di fallback in caso di errore file
        return {{0.0, 0.0}, {2000.0, 100000.0}}; 
    }

    vector<vector<double>> numeri;
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        double val1, val2;
        if (ss >> val1 >> val2) {
            numeri.push_back({val1, val2});
        }
    }
    file.close();
    return numeri;
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
    // Calcolo carichi verticali (Load Transfer Longitudinale)
    // Fz_front = (m*g*lb - m*ax*h) / L
    // Fz_rear  = (m*g*la + m*ax*h) / L
    double L = la + lb;
    double weight_term = m * 9.81;
    double transfer_term = m * acceleration_x * Ycm;

    double Fant = (weight_term * lb - transfer_term) / L;
    double Frear = (weight_term * la + transfer_term) / L;

    // Interpolazione delle stiffness Ka e Kp dalla tabella
    double Ka = interpolate(Fant, numeri);
    double Kp = interpolate(Frear, numeri);

    return make_pair(Ka, Kp);
}