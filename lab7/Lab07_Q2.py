import numpy as np
import matplotlib.pyplot as plt

"""
Exercise 8.13 - Part a
Purpose - Calculate and Plot Earth's orbit

# Pseudo-code
# Set constants
#   'G': Gravitational constant
#   'M': Mass of the Earth
#   'H': Interval = 1week
#   'x_0': Initial x-position
#   'y_0': Initial y-position
#   'vx_0': Initial x-velocity
#   'vy_0': Initial y-velocity
# Define bulirsch_stoer function
#   1) Define modified_midpoint function to compute R_n_1
#   2) Define richardson function to compute R_n_m
# Create empty arrays to hold x, y positions
# Create an array that holds x-position, x-velocity, y-position, y-velocity in order
# Perform for loop
#   1) Store x-position and y-position
#   2) Call bulirsch_stoer function
# Plot
"""
G = 6.6738 * (10**(-11)) * (24 * 365 * 60 * 60) ** 2 # Unit: m^3 kg^-1 s^-2
M = 1.9891 * (10**30)                             # Unit: kg
H = 1/52                                          # 1 week = 1/52 year
x_0 = 1.4710 * (10**11)                           # Unit: m
y_0 = 0
vx_0 = 0
vy_0 = 3.0287 * (10**4) * 60 * 60 * 24 * 365      # Unit: m/year
delta = 1000                                      # Accuracy: 1000m/year
ftsz = 24                                         # Fontsize = 16
def bulirsch_stoer(r, H):

    accuracy = H * delta

    def modified_midpoint(r, n):

        def f(r):

            def fx(x, y):
                radius = np.sqrt(x**2 + y**2)
                return -1 * G * M * x / radius**3

            def fy(x, y):
                radius = np.sqrt(x ** 2 + y ** 2)
                return -1 * G * M * y / radius**3

            x = r[0]
            vx = r[1]
            y = r[2]
            vy = r[3]

            return np.array([vx, fx(x, y), vy, fy(x, y)], float)

        r = np.copy(r)
        h = H/n
        r2 = r + 0.5 * h * f(r)
        r += h * f(r2)
        for i in range(n-1):
            r2 += h * f(r)
            r += h * f(r2)

        return 0.5 * (r + r2 + 0.5 * h * f(r))

    def richardson(R1, n):

        def R_n_m(m):
            return R2[m-2] + (R2[m-2] - R1[m-2])/ ((n / (n-1)) ** (2 * (m-1)) - 1) # Equation 8.103

        # Compute R_n_1
        R2 = [modified_midpoint(r, n)]

        # Compute the rest of the row

        for m in range(2, n+1): # 2 to n
            R2.append(R_n_m(m))

        # Compute the error
        R2 = np.array(R2, float)
        error_x_y = (R2[m-2] - R1[m-2]) / ((n/ (n-1)) ** (2 * (m-1)) - 1) # Equation 8.101
        error_r = np.sqrt(error_x_y[0] ** 2 + error_x_y[2]**2)            # Compute the error of r = sqrt(x^2 + y^2)
        if error_r > accuracy:
            return richardson(R2, n+1)
        else:
            return R2[n-1]

    R_n_1 = np.array([modified_midpoint(r, 1)], float)
    return richardson(R_n_1, 2)

x_earth = []
y_earth = []
r = np.array([x_0, vx_0, y_0, vy_0], float)
for i in range(100):
    x_earth.append(r[0])
    y_earth.append(r[2])
    r = bulirsch_stoer(r, H)

plt.figure(1)
plt.plot(x_earth, y_earth, label="Earth", color="green")
plt.errorbar(0,0,  fmt='o', color='red', label="Sun")
plt.title("The orbit of the Earth - Bulirsch-Stoer Method", fontsize=ftsz)
plt.xlabel("x (in m)", fontsize=ftsz)
plt.ylabel("y (in m)", fontsize=ftsz)
plt.axis('equal')
plt.axvline(0.)
plt.axhline(0.)
plt.tight_layout()
plt.grid()
plt.legend(fontsize=ftsz)

"""
Exercise 8.13 - Part b
Purpose - Calculate and Plot Pluto's orbit

# Pseudo-code
# Set the initial x-position to 4.4368 * 10^12
# Set the initial x-velocity to 6.1218 * 10^3
# Repeat the code in part a
# Plot
"""


x_0 = 4.4368 * (10**12)   # x_0 = 4.4368 * (10**12) / np.sqrt(2) and y_0 = 4.4368 * (10**12) / np.sqrt(2)
#y_0 = 4.4368 * (10**12) / np.sqrt(2)  # Or x_0 = 4.4368 * (10**12) and y_0 = 0
vy_0 = 6.1218 * (10**3) * 60 * 60 * 24 * 365 # Unit: m/year


x_pluto = []
y_pluto = []
r = np.array([x_0, vx_0, y_0, vy_0], float)
H = 1 # Unit: year
delta = 1000 # 1000m/year

for i in range(260):
    x_pluto.append(r[0])
    y_pluto.append(r[2])
    r = bulirsch_stoer(r, H)

# Orbit of Just Pluto
plt.figure(2)
plt.plot(x_pluto, y_pluto, label="Pluto", color="orange")
plt.errorbar(0,0,  fmt='o', color='red', label="Sun")
plt.title("Orbit of Pluto", fontsize=ftsz)
plt.xlabel("x (in m)", fontsize=ftsz)
plt.ylabel("y (in m)", fontsize=ftsz)
#plt.tight_layout()
plt.grid()
plt.axis('equal')
plt.axvline(0.)
plt.axhline(0.)
plt.legend(fontsize=ftsz)


# Orbit of both planets
plt.figure(3)
plt.plot(x_earth, y_earth, label="Earth", color="green")
plt.plot(x_pluto, y_pluto, label="Pluto", color="orange")
plt.errorbar(0,0,  fmt='o', color='red', label="Sun")
plt.title("Orbits of Earth and Pluto", fontsize=ftsz)
plt.xlabel("x (in m)", fontsize=ftsz)
plt.ylabel("y (in m)", fontsize=ftsz)
#plt.tight_layout()
plt.grid()
plt.axis('equal')
plt.axvline(0.)
plt.axhline(0.)
plt.legend(fontsize=ftsz)

plt.show()