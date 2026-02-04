cmake_minimum_required(VERSION 3.10)
project(MPCLongitudinalLateral)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova le librerie necessarie
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED) # Questo include anche OSQP solitamente

# Include directories
include_directories(
    ${EIGEN3_INCLUDE_DIR}
    ${OSQP_EIGEN_INCLUDE_DIRS}
    . 
)

# Crea l'eseguibile di test
# Uniamo tutti i file sorgente
add_executable(test_mpc main.cpp MPC.cpp)

# Linka le librerie all'eseguibile
target_link_libraries(test_mpc
    OsqpEigen::OsqpEigen
    Eigen3::Eigen
)

# Eseguibile per la simulazione
add_executable(run_simulation simulation.cpp MPC.cpp)
target_link_libraries(run_simulation OsqpEigen::OsqpEigen Eigen3::Eigen)