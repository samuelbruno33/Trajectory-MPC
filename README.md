# Coupled Longitudinal & Lateral MPC for EV
UNIFI Master's degree Thesis Work: Extended MPC for Formula Student Autonomous System

## Introduction
This project extends the control architecture of the Firenze Race Team driverless vehicle. It implements a **Coupled Linear Time-Varying Model Predictive Control (LTV-MPC)** algorithm capable of handling both lateral dynamics (steering) and longitudinal dynamics (traction/braking) simultaneously.

Unlike the previous version, which treated velocity as a constant parameter, this controller optimizes the vehicle's speed profile to maximize performance while respecting the **electric powertrain power constraints** ($P_{max} = 80 kW$) and tire friction limits.

## Formulation
The controller uses a dynamic bicycle model formulated in Cartesian coordinates to preserve linearity with respect to control inputs. The model is linearized at each time step around the current operating point.

### State Space
The state vector $x$ consists of 6 variables:
$$x = [X, Y, \theta, v_x, v_y, \omega]^T$$
Where:
* $X, Y$: Global position
* $\theta$: Yaw angle
* $v_x$: Longitudinal velocity (Controlled state)
* $v_y$: Lateral velocity
* $\omega$: Yaw rate

### Control Inputs
The control vector $u$ consists of 2 inputs:
$$u = [F_x, \delta]^T$$
Where:
* $F_x$: Longitudinal Force (Traction positive, Braking negative)
* $\delta$: Steering angle

### Constraints
The optimization problem respects hard constraints on actuators and dynamic constraints on power:
1.  **Steering limits:** $\delta \in [\delta_{min}, \delta_{max}]$
2.  **Actuator Slew Rate:** Limited $\Delta F_x$ and $\Delta \delta$ to ensure smooth control and mechanical integrity.
3.  **Dynamic Power Limit:** The maximum tractive force is bounded by the engine power map:
    $$F_{x,max} = \min(F_{peak}, \frac{P_{max}}{v_x})$$

## Dependencies
This project relies on the following C++ libraries:
* **Eigen3**: For high-performance linear algebra operations.
* **OSQP**: For solving the quadratic programming (QP) optimization problem.
* **osqp-eigen**: A C++ wrapper to interface Eigen with OSQP.

### Installation (Windows Testing)
Dependencies are managed via **vcpkg**. To set up the environment:

1.  Clone vcpkg (if not already present):
    ```powershell
    git clone [https://github.com/microsoft/vcpkg.git](https://github.com/microsoft/vcpkg.git)
    .\vcpkg\bootstrap-vcpkg.bat
    ```
2.  Install libraries:
    ```powershell
    .\vcpkg\vcpkg install eigen3 osqp osqp-eigen
    ```

## Build & Run
To compile the standalone test on Windows using CMake:

```powershell
mkdir build
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=[PATH_TO_VCPKG]/scripts/buildsystems/vcpkg.cmake
cmake --build .
```

Run the test executable:
```powershell
.\Debug\test_mpc.exe
```

## Credits
Based on the Bachelor's Thesis work by Federico Monetti (Lateral Control). Extended by Samuel Bruno (Longitudinal & Lateral Control for EV).