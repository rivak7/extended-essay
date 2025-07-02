import math
import csv

g = 9.80065                # Gravity (m/s^2)
Cd = 0.67                  # Drag coefficient (dimensionless)
A = 0.003425               # Cross-sectional area (m^2)
p = 1.229                  # Air density (kg/m^3)
T = 73.4                   # Thrust (N)
mB = 0.625                 # Mass during burn phase (kg)
mC = 0.610                 # Mass during coast phase (kg)
tB = 0.85                  # Burn time (s)

a0 = (T) / (mB * g) - 1
b0 = (mB * g) / (0.5 * Cd * A * p)
aC = -1
bC = (mC * g) / (0.5 * Cd * A * p)

vB = 0
a = 0
v = 0
s = 0
t = 0
dt = 0.001

time = []
altitude = []
velocity = []
acceleration = []

# Burn phase
while t <= tB:
    v = math.sqrt(b0 * a0) * math.tanh(g * t * math.sqrt(a0 / b0))
    s = (b0 / g) * math.log(math.cosh(g * t * math.sqrt(a0 / b0)))
    a = T / mB - g - (0.5 * Cd * A * p * v**2 / mB)

    time.append(t)
    altitude.append(s)
    velocity.append(v)
    acceleration.append(a)

    vB = v
    sB = s
    t += dt

# Coast phase
tC = (math.sqrt(bC) / g) * math.atan(vB / math.sqrt(bC)) + tB
sC = (bC / (2 * g)) * math.log(1 + vB**2 / bC)

while t <= tC:
    v = math.sqrt(bC) * math.tan(math.atan(vB / math.sqrt(bC)) - (g * (t - tB)) / math.sqrt(bC))
    s = sC + (bC / g) * math.log(math.cos((g / math.sqrt(bC)) * (tC - t))) + sB
    a = -g - (0.5 * Cd * A * p * v**2 / mC)

    time.append(t)
    altitude.append(s)
    velocity.append(v)
    acceleration.append(a)
    t += dt

with open("analytical.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Time (s)", "Altitude (m)", "Velocity (m/s)", "Acceleration (m/s^2)"])
    for i in range(len(time)):
        writer.writerow([time[i], altitude[i], velocity[i], acceleration[i]])
