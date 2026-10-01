# Model Rocket Apogee Prediction: Analytical vs. Numerical Methods

Code and data for my IB Extended Essay in Mathematics (2024–25), "Comparing the accuracy of analytical and numerical methods for predicting a model rocket's apogee using 1-dimensional kinematics." I wrote it alongside competing in the 2024 American Rocketry Challenge (ARC).

The essay derives the closed form Fehskens-Malewicki solution for a rocket's altitude during the burn and coast phases, then compares it against two numerical integrators: explicit Euler and fourth-order Runge–Kutta (RK4), using onboard altimeter data from a real flight. The PDF is included exactly as submitted. See the [retrospective](#retrospective-2026) below for problems I've since found and corrected.

Flight data courtesy of Nishant Vikramaditya.

## Contents

| File | Description |
|---|---|
| `Extended_Essay.pdf` | The essay, as submitted |
| `analytical_implement.py` | Analytical solution with constant thrust, mass, and air density in each (burn and coast) phase |
| `numerical_implement.py` | Explicit Euler and RK4 integration with a time-varying thrust curve, mass, and air density |
| `interpolate.py` | Resamples flight data to a 1 ms timestep |
| `AeroTech_F67W.csv` | Motor thrust curve, from [thrustcurve.org](https://www.thrustcurve.org/motors/AeroTech/F67W/) |
| `Data.xlsx` | Flight data and each method's predicted altitude, velocity, and acceleration |
| `Method_Comparison_Graph.png` | Figure 4 from the essay |

## Running

The scripts use only the Python standard library.

```bash
python analytical_implement.py      # writes analytical.csv
```

For the numerical methods, uncomment a line in the driver at the bottom of `numerical_implement.py`, then run it. Each solver writes `rocket_trajectory.csv` and returns the predicted apogee.

The raw altimeter export that `interpolate.py` reads (`Real.csv`) is not included. The interpolated flight data is in the "Experimental" columns of `Data.xlsx`.

## Results as submitted

| Method | Apogee (m) |
|---|---|
| Flight data | 254.022 |
| Analytical | 257.365 |
| Explicit Euler | 260.542 |
| RK4 | 257.099 |

The essay concluded that RK4 was the most accurate method and the analytical method was surprisingly competitive.

## Retrospective (2026)

Rereading this code in 2026, I found a bug that produced the essay's main result, along with a major flaw in experimental design...

### RK4 implementation error

The original `update_rk4` computed intermediate states for each stage (`s_temp`, `v_temp`, `t_temp`), but `acceleration()` and the air-density update read the *current* state instead. So stages k_2 through k_4 were all evaluated at the start of the step rather than at their intermediate points. Separately, the rocket's mass was reset to its starting value at the end of every step, so propellant never burned: the rocket stayed at its 0.640 kg wet mass for the whole flight.

The fixed version computes derivatives as a pure function of the state (altitude, velocity, mass), so each stage uses its own intermediate state and mass is integrated along with altitude and velocity.

### Corrected results

| Method | Apogee (m) | Error (m) |
|---|---|---|
| Flight data | 254.022 | -- |
| Analytical | 257.264 | +3.24 |
| Explicit Euler | 260.522 | +6.50 |
| RK4 | 260.358 | +6.34 |

all at dt = 1 ms

The submitted "RK4" landed closest to the flight data only because the bug kept the rocket about 30 g too heavy, which happened to offset the overprediction from ignoring wind. Implemented correctly, RK4 and explicit Euler agree to within 0.2 m.

### Experimental design

Comparing each method against a single real flight mostly measures **modeling error**: how well the physical model (an estimated drag coefficient of 0.67, no wind, 1D motion) matches reality. It says little about **numerical error**, which at dt = 1 ms is a small fraction of a meter for both integrators:

| dt (s) | Explicit Euler (m) | RK4 (m) |
|---|---|---|
| 0.01 | 261.974 | 260.363 |
| 0.001 | 260.522 | 260.358 |
| 0.0001 | 260.374 | 260.358 |

Euler's error shrinks roughly tenfold for each tenfold decrease in dt, as expected for a first-order method, while RK4 has already converged at dt = 0.01 s. Apogee here is the highest sampled altitude, so it is resolved only to within one timestep.

The differences between methods are much smaller than the uncertainty in the model's inputs, so one flight can't rank them. A better design would separate the two error sources:

- **Numerical error:** run a convergence study like the one above, and compare the integrators against the analytical solution under identical constant-parameter assumptions, where the analytical answer is exact.
- **Modeling error:** test sensitivity to uncertain inputs like the drag coefficient, and compare against multiple flights.

### Other fixes

- **Explicit Euler:** The original Euler update used the *new* velocity to update altitude, which is semi-implicit Euler, not the explicit method described in equation (18). It now uses only the state at t_n. This changes the 1 ms apogee from 260.542 m to 260.522 m.
- **Gravity constant:** `analytical_implement.py` used g = 9.80065 instead of 9.80665 m/s^2. Fixing it changes the analytical apogee from 257.365 m to 257.264 m.
- **Clarity:** `J` was commented as specific impulse but is the motor's total impulse (N*s), which is what the mass-flow formula dm/dt = −T\*m_prop/J needs. The hardcoded propellant mass is now a parameter, `m_prop`.

The code exactly as submitted is preserved at the [`ee-submitted`](https://github.com/rivak7/extended-essay/tree/ee-submitted) tag.
