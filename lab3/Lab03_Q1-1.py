import numpy as np
from scipy import special
import matplotlib.pyplot as plt
from gaussxw import gaussxw

"""
Q1 Part A.i
Purpose - Evaluate the Dawson function using Trapezoidal Rule, Simpson's rule, and scipy.special.dawsn

# Pseudocode
# Import scipy.special module
# Define the dawson function
# Create the starting point and end point (a, b)
# Define the number of slices (N)
# Calculate the width of each slice (h = (b-a)/N)
# Use equation 5.3 to compute the Dawson integral using Trapezoidal Rule
# Use equation 5.9 to compute the Dawson integral using Simpson's Rule
# Use scipy.special.dawsn to compute the Dawson integral
# Output the value of the evaludated integral with the number of slices and the method used
"""


def f(t):
    return np.exp(t**2)

N = np.linspace(8, 2048, 100)                   # 8 slices
N = N.astype(int)
a = 0.0                 # Start point
b = 4.0                 # End point
h = (b-a)/N             # Width of a slice

# Trapezoidal Rule

trapezoidal = np.zeros(len(N))

for i in range(len(N)):
    s_trapezoidal = 0.5 * f(a) + 0.5 * f(b)
    for k in range(1, N[i]):
        s_trapezoidal += f(a + k * h[i])

    trapezoidal[i] = np.exp((-1)*(b**2)) * h[i] * s_trapezoidal

    print("Part A.i - Trapezoidal rule - ", N[i], " number of slices : ", trapezoidal[i]) # Equation 5.3

# Simpson's Rule

simpson = np.zeros(len(N))
for i in range(len(N)):
    s_simpson = f(a) + f(b)

    for k in range(1, N[i], 2):
        s_simpson += 4 * f(a + k * h[i])

    for k in range(2, N[i], 2):
        s_simpson += 2 * f(a + k * h[i])

    simpson[i] = np.exp((-1) * (b ** 2)) * (1 / 3) * h[i] * s_simpson
    print("Part A.i - Simpson's rule - ", N[i], " number of slices : ", simpson[i])  # Equation 5.9

# Scipy.special.dawsn
print("Part A.i - Scipy.special.dawsn : ", special.dawsn(b))

print("===============================================================================")
"""
Q1 Part A.ii
Purpose - Calculate the relative error compared to the true value of D(4)

# Pseudocode
# Define a new number of slices by multiplying the original one by two
# Calculate the width of new slices
# Repeat Part A.i with the new slices and width
# Compute the integral using trapezoidal rule and simpson's rule
# Subtract the original value from the new value = relative error
"""
N2 = 2*N
h2 = (b-a)/N2

# Relative Error - Trapezoidal Rule

trapezoidal2 = np.zeros(len(N2))
error_trapezoidal = np.zeros(len(N2))

for i in range(len(N2)):
    s_trapezoidal = 0.5 * f(a) + 0.5 * f(b)
    for k in range(1, N2[i]):
        s_trapezoidal += f(a + k * h2[i])

    trapezoidal2[i] = np.exp((-1) * (b ** 2)) * h2[i] * s_trapezoidal
    error_trapezoidal[i] = trapezoidal2[i] - trapezoidal[i]

    print("Trapezoidal - Relative error at N = ", N[i], " : ", error_trapezoidal[i])

# Relative Error - Simpsons's Rule

simpson2 = np.zeros(len(N2))
error_simpson = np.zeros(len(N2))

for i in range(len(N2)):
    s_simpson = f(a) + f(b)

    for k in range(1, N2[i], 2):
        s_simpson += 4 * f(a + k * h2[i])

    for k in range(2, N2[i], 2):
        s_simpson += 2 * f(a + k * h2[i])

    simpson2[i] = np.exp((-1) * (b ** 2)) * (1 / 3) * h2[i] * s_simpson
    error_simpson[i] = simpson2[i] - simpson[i]

    print("Simpson - Relative error at N = ", N[i], " : ", error_simpson[i])

# Plot

plt.figure(1)
plt.loglog(N, np.abs(error_trapezoidal), '.', label="Relative Error - Trapezoidal")
plt.loglog(N, error_simpson, '.', label="Relative Error - Simpson")
plt.xlabel("log(Number of slices)", fontsize=16)
plt.ylabel("log(Relative Error)", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title("Number of slices vs Relative Error in log-log plot", fontsize=16)
plt.legend(fontsize=16)



"""
Q1 Part B
Purpose - Calculate the probability of blowing snow using a Gaussian quadrature with N=100

# Pseudocode
# Define the number of slices
# Define the starting and ending points
# Define the initial conditions - wind strength and the age of the snow
# Create an array of averge hourly temperature
# Define the function to calculate equation2
# Apply the change of variables formula
# Perform gaussian quadrature 
# Plot
"""

N = 100
a = 0
b = np.array([6, 8, 10]) # Average hourly windspeed at a height of 10m, Unit : m/s^-1
th = np.array([24, 48, 72]) # Snow surface age, Unit: hours

ta = np.linspace(-100, 100, 100) # # Average hourly temperature
ta = ta.astype(int)

def delta(ta):              # Standard Deviation
    return 4.3+0.145*ta+0.00196*(ta**2)

def g0(u): # th = 24
    return np.exp((-1)*((11.2+0.365*ta+0.00706*(ta**2)+0.9*np.log(th[0])-u)**2)/(2*(delta(ta))**2))

def g1(u): # th = 48
    return np.exp((-1)*((11.2+0.365*ta+0.00706*(ta**2)+0.9*np.log(th[1])-u)**2)/(2*(delta(ta))**2))

def g2(u): # th = 72
    return np.exp((-1)*((11.2+0.365*ta+0.00706*(ta**2)+0.9*np.log(th[2])-u)**2)/(2*(delta(ta))**2))

x, w = gaussxw(N)
xp0 = 0.5*(b[0]-a)*x + 0.5*(b[0]+a)     # u10 = 6
xp1 = 0.5*(b[1]-a)*x + 0.5*(b[1]+a)     # u10 = 8
xp2 = 0.5*(b[2]-a)*x + 0.5*(b[2]+a)     # u10 = 10

wp0 = 0.5*(b[0]-a)*w    # u10 = 6
wp1 = 0.5*(b[1]-a)*w    # u10 = 8
wp2 = 0.5*(b[2]-a)*w    # u10 = 10

#  Integration
# u10 = 6, th = 24
s0 = 0.0
for k in range(N):
    s0 += wp0[k]*g0(xp0[k])

# u10 = 6, th = 48
s1 = 0.0
for k in range(N):
    s1 += wp0[k]*g1(xp0[k])

# u10 = 6, th = 72
s2 = 0.0
for k in range(N):
    s2 += wp0[k]*g2(xp0[k])

# u10 = 8, th = 24
s3 = 0.0
for k in range(N):
    s3 += wp1[k]*g0(xp1[k])

# u10 = 8, th = 48
s4 = 0.0
for k in range(N):
    s4 += wp1[k]*g1(xp1[k])

# u10 = 8, th = 72
s5 = 0.0
for k in range(N):
    s5 += wp1[k]*g2(xp1[k])

# u10 = 10, th = 24
s6 = 0.0
for k in range(N):
    s6 += wp2[k]*g0(xp2[k])

# u10 = 10, th = 48
s7 = 0.0
for k in range(N):
    s7 += wp2[k]*g1(xp2[k])

# u10 = 10, th = 72
s8 = 0.0
for k in range(N):
    s8 += wp2[k]*g2(xp2[k])

# Multiply Inverse of sqrt(2pi)*standard deviation
for i in range(len(ta)):
    s0[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s0[i]
    s1[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s1[i]
    s2[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s2[i]
    s3[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s3[i]
    s4[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s4[i]
    s5[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s5[i]
    s6[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s6[i]
    s7[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s7[i]
    s8[i] = (1/np.sqrt(2*np.pi)) * (1/delta(ta[i])) * s8[i]

plt.figure(3)
plt.plot(ta, s6, label="u$_{10}$ = 10 m/s, t$_{h}$ = 24 hrs")
plt.plot(ta, s7, label="u$_{10}$ = 10 m/s, t$_{h}$ = 48 hrs")
plt.plot(ta, s8, label="u$_{10}$ = 10 m/s, t$_{h}$ = 72 hrs")
plt.plot(ta, s3, label="u$_{10}$ = 8 m/s, t$_{h}$ = 24 hrs")
plt.plot(ta, s4, label="u$_{10}$ = 8 m/s, t$_{h}$ = 48 hrs")
plt.plot(ta, s5, label="u$_{10}$ = 8 m/s, t$_{h}$ = 72 hrs")
plt.plot(ta, s0, label="u$_{10}$ = 6 m/s, t$_{h}$ = 24 hrs")
plt.plot(ta, s1, label="u$_{10}$ = 6 m/s, t$_{h}$ = 48 hrs")
plt.plot(ta, s2, label="u$_{10}$ = 6 m/s, t$_{h}$ = 72 hrs")
plt.xlabel("T$_{a}$ - Average Hourly Temperature (in $^\circ$C)", fontsize=16)
plt.ylabel("P(u$_{10}$, T$_{a}$, t$_{h}$) - Probability of blowing snow", fontsize=16)
plt.title("Probability of blowing snow vs Average Hourly Tempearture", fontsize=16)
plt.legend(fontsize=16)



plt.show()