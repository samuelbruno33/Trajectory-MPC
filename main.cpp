#include <iostream>
#include <vector>
#include "MPC.hpp"
#include "geometry.hpp"
#include <Eigen/Dense>

int main() {
    // 1. Istanzia il controller
    MPC mpc_controller;
    
    // 2. Crea uno stato finto iniziale
    // [X, Y, Theta, Vx, Vy, Omega]
    Eigen::VectorXd x0(6);
    x0 << 0.0,   // X
          0.0,   // Y (Siamo al centro della pista)
          0.0,   // Theta (Allineati)
          10.0,  // Vx = 10 m/s (36 km/h)
          0.0,   // Vy
          0.0;   // Omega

    // 3. Crea una traiettoria di riferimento (Una retta lungo l'asse X)
    std::vector<Point> waypoints;
    for (int i = 0; i < 30; ++i) {
        // Waypoints ogni 1 metro lungo X, con Y=0
        waypoints.push_back(Point(static_cast<double>(i), 0.0));
    }

    // 4. Aggiorna il modello (Linearizzazione & Discretizzazione)
    // Supponiamo acc_y misurata = 0
    mpc_controller.updateDiscretization(x0, 0.0);

    // 5. Calcola il controllo
    // Reference speed = 15 m/s (vogliamo accelerare)
    double v_ref = 15.0; 
    std::pair<double, double> result = mpc_controller.compute(x0, waypoints, v_ref);

    // 6. Stampa i risultati
    std::cout << "--- RISULTATI TEST MPC ---" << std::endl;
    std::cout << "Stato Iniziale Vx: " << x0(3) << " m/s" << std::endl;
    std::cout << "Target Vx: " << v_ref << " m/s" << std::endl;
    std::cout << "Input Calcolati:" << std::endl;
    std::cout << ">> Forza Longitudinale (Fx): " << result.first << " N" << std::endl;
    std::cout << ">> Angolo Sterzo (Delta):    " << result.second << " rad" << std::endl;

    // Check rapido
    if(result.first > 0) std::cout << "[OK] Il veicolo sta accelerando come previsto." << std::endl;
    else std::cout << "[?] Il veicolo frena o è in coasting." << std::endl;

    return 0;
}