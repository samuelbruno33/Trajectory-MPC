# Trajectory MPC for EV
UNIFI Bachelor's degree Thesis Work: Extended MPC for Formula Student Autonomous System

## Introduction
This project extends the control architecture of the Firenze Race Team driverless vehicle. It implements a **Linear Time-Varying Model Predictive Control (LTV-MPC)** algorithm capable of handling both lateral dynamics (steering) and longitudinal dynamics (traction/braking) simultaneously.

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


## Installation & Dependencies
The project uses **vcpkg** for dependency management on Windows.

1.  **Clone vcpkg:**
    ```powershell
    git clone [https://github.com/microsoft/vcpkg.git](https://github.com/microsoft/vcpkg.git)
    .\vcpkg\bootstrap-vcpkg.bat
    ```
2.  **Install Libraries:**
    ```powershell
    .\vcpkg\vcpkg install eigen3 osqp osqp-eigen
    ```

## How to Run the Simulation
Follow these steps to reproduce the thesis results on Windows:

### Step A: Build
```powershell
mkdir build
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=[PATH_TO_VCPKG]/scripts/buildsystems/vcpkg.cmake
cmake --build .

```

### Step B: Run Simulation

Ensure the file `cornering_stiffness_vs_vertical_load.txt` is present in the executable directory.

```powershell
.\Debug\run_simulation.exe

```

*This generates a `comparison_results.csv` file containing the telemetry of both controllers.*

### Step C: Generate Plots

Use the provided Python script to visualize the comparison:

```powershell
python conversion_thesis.py

```

*This generates `grafico_1_traiettoria.png`, `grafico_2_longitudinale.png` and `grafico_3_dinamica_longitudinale.png`.*

## Experimental Results

The following statistics were collected during a comparative test on the "Chicane" track scenario.

### Performance Metrics (RMSE)

The Coupled MPC demonstrates superior tracking capabilities compared to the PID+PP Baseline.

| Metric | MPC (Coupled) | PID (Baseline) | Improvement |
| --- | --- | --- | --- |
| **Lateral Error (RMSE)** | **1.381 m** | 6.213 m | **+77.8%** 🟢 |
| **Velocity Error (RMSE)** | **7.232 m/s** | 16.513 m/s | **+56.2%** 🟢 |

### Power & Energy Analysis

* **Energy Consumed:**
* MPC: 187.8 kJ
* PID: 153.7 kJ
* *Note: The PID consumed less energy simply because it failed to reach the target speeds (high velocity error), whereas the MPC maximized performance.*


* **Constraint Satisfaction:**
* The MPC successfully managed the power limitation, saturating the actuator only when necessary to respect the  limit curve.
* The PID baseline showed inability to anticipate the curve, resulting in massive tracking errors.


## Credits

* **Federico Monetti:** Original Lateral MPC implementation.
* **Samuel Bruno:** Coupled Longitudinal/Lateral formulation and EV Constraints.