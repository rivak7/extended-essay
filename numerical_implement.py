import csv
import math

class Rocket:
    def __init__(self, thrust_file, Cd, m, A, J, psea, p0, T0, g, output_file):
        self.s = 0.01  # Altitude (m)
        self.v = 0  # Velocity (m/s)
        self.t = 0  # Time (s)
        self.dt = 0  # Time step (s)
        self.m = m  # Initial mass (kg)
        self.A = A  # Cross-sectional area (m^2)
        self.Cd = Cd  # Drag coefficient
        self.J = J  # Specific impulse (s)
        self.g = g  # Gravity (m/s^2)
        self.psea = psea  # Sea level air density (kg/m^3)
        self.p0 = p0  # Sea level pressure (Pa)
        self.T0 = T0  # Sea level temperature (K)
        self.p = psea  # Current air density (kg/m^3)
        self.thrustcurve = self.load_thrust_curve(thrust_file)
        self.output_file = output_file  # CSV output file

    def load_thrust_curve(self, file):
        thrustcurve = []
        with open(file, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                try:
                    thrustcurve.append((float(row[0]), float(row[1])))
                except ValueError:
                    continue  # Skip headers or invalid rows
        return thrustcurve

    def get_thrust(self, t):
        for i in range(len(self.thrustcurve)):
            if t < self.thrustcurve[i][0]:
                if i == 0:
                    return self.thrustcurve[i][1]
                t1, T1 = self.thrustcurve[i - 1]
                t2, T2 = self.thrustcurve[i]
                return ((t - t1) / (t2 - t1)) * (T2 - T1) + T1
        return 0  # No thrust after the last point

    def update_pressure(self):
        self.p = ((0.0289652 * self.p0) / (8.31446 * self.T0)) * math.pow((1 - (0.0065 * self.s / self.T0)),
                                                                           (((0.0289652 * self.g) / (8.31446 * 0.0065)) - 1))

    def update_mass(self):
        if self.get_thrust(self.t)>=0:
            dm = (self.get_thrust(self.t) * self.dt * 0.030) / self.J
            self.m -= dm

    def acceleration(self):
        drag = 0.5 * self.Cd * self.A * self.p * self.v**2 / self.m
        thrust = self.get_thrust(self.t) / self.m
        return thrust - self.g - drag

    def update_euler(self):
        self.update_pressure()
        self.update_mass()
        a = self.acceleration()
        self.v += a * self.dt
        self.s += self.v * self.dt
        self.t += self.dt

        # Output the current state to CSV
        with open(self.output_file, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([self.t, self.s, self.v, a])

    def solve_euler(self, step):
        self.dt = step

        # Write headers to the CSV file
        with open(self.output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Time (s)", "Altitude (m)", "Velocity (m/s)", "Acceleration (m/s^2)"])

        while self.v >= 0:
            self.update_euler()

    def update_rk4(self):
        # Save the initial state
        initial_mass = self.m
        initial_pressure = self.p

        # Update pressure and mass at the start of the timestep
        self.update_pressure()
        self.update_mass()

        # Calculate k1
        k1v = self.acceleration()
        k1s = self.v

        # Intermediate step for k2
        self.t_temp = self.t + 0.5 * self.dt  # Half-step for time
        self.s_temp = self.s + 0.5 * k1s * self.dt  # Half-step for position
        self.v_temp = self.v + 0.5 * k1v * self.dt  # Half-step for velocity
        self.update_pressure()
        self.update_mass()
        k2v = self.acceleration()
        k2s = self.v

        # Intermediate step for k3
        self.s_temp = self.s + 0.5 * k2s * self.dt - 0.5 * k1s * self.dt  # Correct for previous half-step
        self.v_temp = self.v + 0.5 * k2v * self.dt - 0.5 * k1v * self.dt  # Correct for previous half-step
        self.update_pressure()
        self.update_mass()
        k3v = self.acceleration()
        k3s = self.v

        # Final step for k4
        self.s_temp = self.s + k3s * self.dt - 0.5 * k2s * self.dt  # Correct for previous half-step
        self.v_temp = self.v + k3v * self.dt - 0.5 * k2v * self.dt  # Correct for previous half-step
        self.update_pressure()
        self.update_mass()
        k4v = self.acceleration()
        k4s = self.v

        # Update state variables using weighted average
        self.v += (self.dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
        self.s += (self.dt / 6) * (k1s + 2 * k2s + 2 * k3s + k4s)

        # Restore mass and pressure to the correct state
        self.m = initial_mass
        self.p = initial_pressure

        # Advance time by the full timestep
        self.t += self.dt  # Now we update time once at the end

        # Write to CSV every timestep (0.001 s)
        with open(self.output_file, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([self.t, self.s, self.v, self.acceleration()])

    def solve_rk4(self, step):
        self.dt = step

        # Write headers to the CSV file
        with open(self.output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Time (s)", "Altitude (m)", "Velocity (m/s)", "Acceleration (m/s^2)"])

        while self.v >= 0:
            self.update_rk4()

if __name__ == "__main__":
    rocket = Rocket(
        thrust_file="AeroTech_F67W.csv",
        Cd=0.67,
        m=0.640,
        A=0.003425,
        J=61.1,
        psea=1.229,
        p0=101625,
        T0=296,
        g=9.80665,
        output_file="rocket_trajectory.csv"  # Specify the output CSV file
    )

    # Driver code: call either function to implement either explicit Euler or RK4
    # rocket.solve_euler(0.001)
	# rocket.solve_RK4(0.001)
	