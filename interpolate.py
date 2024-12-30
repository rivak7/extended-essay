import csv

dt = 0.001

def load_data(file, column):
    data = []
    with open(file, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            try:
                data.append((float(row[0]), float(row[column])))
            except ValueError:
                continue  # Skip headers
    return data

s_sensor = load_data("Real.csv", 1)
v_sensor = load_data("Real.csv", 2)
a_sensor = load_data("Real.csv", 3)

# Get the maximum time from the first column of the sensor data
t_max = max(s_sensor, key=lambda x: x[0])[0]

def interpolate_data(column, dt, t_max):
    interpolated = []
    t = 0
    i = 0
    while t <= t_max:
        # Find the two points to interpolate between
        while i < len(column) - 1 and column[i + 1][0] < t:
            i += 1
        if i == len(column) - 1:
            # If we are at the end of the dataset, append the last value
            interpolated.append(column[i][1])
        else:
            # Interpolate between column[i] and column[i + 1]
            t1, T1 = column[i]
            t2, T2 = column[i + 1]
            interpolated.append(((t - t1) / (t2 - t1)) * (T2 - T1) + T1)
        t += dt
    return interpolated

s = interpolate_data(s_sensor, dt, t_max)
v = interpolate_data(v_sensor, dt, t_max)
a = interpolate_data(a_sensor, dt, t_max)

def save_interpolated_data(file, t_max, dt, s, v, a):
    with open(file, "w", newline="") as f:
        writer = csv.writer(f)
        # Write the header row
        writer.writerow(["Time", "s", "v", "a"])
        # Write the data rows
        t = 0
        for i in range(len(s)):
            writer.writerow([round(t, 3), s[i], v[i], a[i]])
            t += dt

# Save to a CSV file
save_interpolated_data("interpolateddata.csv", t_max, dt, s, v, a)