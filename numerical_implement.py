import csv
import math

M_AIR = 0.0289652      # Molar mass of dry air (kg/mol)
R_GAS = 8.31446        # Universal gas constant (J/(mol*K))
LAPSE = 0.0065         # Temperature lapse rate (K/m)


class Rocket:
    def __init__(self, thrust_file, Cd, m, A, J, m_prop, p0, T0, g, output_file):
        self.s0 = 0.01            # Initial altitude (m)
        self.v0 = 0.0             # Initial velocity (m/s)
        self.m0 = m               # Initial (wet) mass (kg)
        self.A = A                # Cross-sectional area (m^2)
        self.Cd = Cd              # Drag coefficient
        self.J = J                # Motor total impulse (N*s)
        self.m_prop = m_prop      # Propellant mass (kg)
        self.g = g                # Gravity (m/s^2)
        self.p0 = p0              # Launch site pressure (Pa)
        self.T0 = T0              # Launch site temperature (K)
        self.thrustcurve = self.load_thrust_curve(thrust_file)
        self.output_file = output_file

    def load_thrust_curve(self, file):
        thrustcurve = []
        with open(file, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                try:
                    thrustcurve.append((float(row[0]), float(row[1])))
                except ValueError:
                    continue
        return thrustcurve

    def get_thrust(self, t):
        """Linearly interpolate the thrust curve; zero after burnout."""
        for i in range(len(self.thrustcurve)):
            if t < self.thrustcurve[i][0]:
                if i == 0:
                    return self.thrustcurve[i][1]
                t1, T1 = self.thrustcurve[i - 1]
                t2, T2 = self.thrustcurve[i]
                return ((t - t1) / (t2 - t1)) * (T2 - T1) + T1
        return 0.0

    def air_density(self, s):
        """Barometric formula for air density (kg/m^3) at altitude s (m)."""
        exponent = (M_AIR * self.g) / (R_GAS * LAPSE) - 1
        return ((M_AIR * self.p0) / (R_GAS * self.T0)) * (1 - LAPSE * s / self.T0) ** exponent

    def derivatives(self, t, s, v, m):
        """
        RHS of the ODE system for the state (s, v, m).
        Every quantity is computed from the arguments so each RK4 stage uses
        its own intermediate state.
        """
        thrust = self.get_thrust(t)
        drag = 0.5 * self.Cd * self.A * self.air_density(s) * v**2
        ds_dt = v
        dv_dt = (thrust - m * self.g - drag) / m
        # Mass flow is proportional to thrust: dm/dt = -T * m_prop / J
        dm_dt = -thrust * self.m_prop / self.J
        return ds_dt, dv_dt, dm_dt

    def _solve(self, step_fn, dt):
        t, s, v, m = 0.0, self.s0, self.v0, self.m0
        rows = []
        while v >= 0:
            s, v, m = step_fn(t, s, v, m, dt)
            t += dt
            rows.append([t, s, v, self.derivatives(t, s, v, m)[1]])

        with open(self.output_file, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Time (s)", "Altitude (m)", "Velocity (m/s)", "Acceleration (m/s^2)"])
            writer.writerows(rows)
        return max(row[1] for row in rows)  # apogee (m)

    def euler_step(self, t, s, v, m, dt):
        """Explicit Euler: every update uses only the state at t_n (eqs 17-18)"""
        ds, dv, dm = self.derivatives(t, s, v, m)
        return s + ds * dt, v + dv * dt, m + dm * dt

    def rk4_step(self, t, s, v, m, dt):
        """Classical fourth-order Runge-Kutta on the full state (s, v, m) (eq 19)"""
        k1 = self.derivatives(t, s, v, m)
        k2 = self.derivatives(t + dt / 2, s + dt / 2 * k1[0], v + dt / 2 * k1[1], m + dt / 2 * k1[2])
        k3 = self.derivatives(t + dt / 2, s + dt / 2 * k2[0], v + dt / 2 * k2[1], m + dt / 2 * k2[2])
        k4 = self.derivatives(t + dt, s + dt * k3[0], v + dt * k3[1], m + dt * k3[2])
        return tuple(
            y + dt / 6 * (a + 2 * b + 2 * c + d)
            for y, a, b, c, d in zip((s, v, m), k1, k2, k3, k4)
        )

    def solve_euler(self, dt):
        return self._solve(self.euler_step, dt)

    def solve_rk4(self, dt):
        return self._solve(self.rk4_step, dt)


# Driver code: call either function to implement explicit Euler or RK4
if __name__ == "__main__":
    rocket = Rocket(
        thrust_file="AeroTech_F67W.csv",
        Cd=0.67,
        m=0.640,
        A=0.003425,
        J=61.1,
        m_prop=0.030,
        p0=101625,
        T0=296,
        g=9.80665,
        output_file="rocket_trajectory.csv",
    )

    # print(f"Explicit Euler apogee: {rocket.solve_euler(0.001):.3f} m")
    # print(f"RK4 apogee: {rocket.solve_rk4(0.001):.3f} m")
