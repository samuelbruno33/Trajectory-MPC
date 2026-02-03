#ifndef MPC_HPP
#define MPC_HPP

#define _USE_MATH_DEFINES

#include <vector>
#include <fstream>
#include <tuple>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <OsqpEigen/OsqpEigen.h>
#include "geometry.hpp"

using namespace std;

// Indici per leggibilità stato e input
enum StateIdx { IDX_X=0, IDX_Y, IDX_THETA, IDX_VX, IDX_VY, IDX_OMEGA };
enum InputIdx { IDX_FX=0, IDX_DELTA };

class MPC {
    public:
        MPC(); 

        // Aggiorna matrici A, B linearizzate attorno allo stato attuale
        // vx è ora parte dello stato, ma lo passiamo per comodità di linearizzazione
        void updateDiscretization(const Eigen::VectorXd& current_state, double acc_y_measured);

        // Calcola il controllo ottimo. Restituisce coppia {Fx, Delta}
        pair<double, double> compute(const Eigen::VectorXd& x0, const vector<Point>& waypoints, double v_ref);

    private:
        vector<vector<double>> load_data();
        pair<double,double> load_transfer(double acc_x, const vector<vector<double>> &numeri);

        // Matrici del modello
        Eigen::MatrixXd A, B;   // Continue
        Eigen::MatrixXd Ad, Bd; // Discrete 
        
        // Parametri Veicolo
        double la = 0.792; 
        double lb = 0.758; 
        double Iz = 98.03; 
        double m = 320; 
        double Ycm = 0.362;
        
        // Parametri Powertrain EV
        double P_max = 80000.0; // 80 kW
        double F_peak = 1500.0; // 1500 N (limite grip/motore a bassa velocità)

        // Dimensioni Problema
        int nx = 6; // [X, Y, Theta, Vx, Vy, Omega]
        int nu = 2; // [Fx, Delta]
        int op = 20; 
        int oc = 10; 
        double dt = 0.05; // Aumentato leggermente per orizzonte più lungo (1s)

        // Pesi Costo
        Eigen::MatrixXd Q;       // Pesi stato
        Eigen::MatrixXd R;       // Pesi input
        Eigen::MatrixXd R_delta; // Pesi slew rate

        // Limiti Costanti
        Eigen::VectorXd xmin, xmax; 
        Eigen::VectorXd umin, umax; // umax sarà aggiornato dinamicamente per Fx
        Eigen::VectorXd u_prev;     // Warm start: [Fx_prev, Delta_prev, ...]

        OsqpEigen::Solver solver;
        bool solver_initialized;

        std::vector<std::vector<double>> numeri; // Tabella cornering stiffness
};

#endif