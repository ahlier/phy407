import numpy as np
from scipy import special
import time
import matplotlib.pyplot as plt

def f(t):
    return np.exp(t**2)

"""
Q2 Part A.i
Purpose - Evaluate the Dawson function at x = 4 and N = 8 slices, using trapezoidal and Simpson's rules and compare
          the result to exact value given by scipy.special.dawsn

# Pseudocode
# Define the number of slices and start & ending points
# Calculate the width of the slice
# Employ equation 5.3 with for loop
"""
N = 8                   # 8 slices
a = 0.0                 # Start point
b = 4.0                 # End point
h = (b-a)/N             # Width of a slice

# Trapezoidal rule

s_trapezoidal = 0.5* f(a) + 0.5 * f(b)
for k in range(1, N):
    s_trapezoidal += f(a+k*h)

print("Part A.i - Trapezoidal rule : ", np.exp((-1)*(b**2)) * h * s_trapezoidal) # Equation 5.3

# Simpson's Rule

s_simpson = f(a) + f(b)

for k in range(1, N, 2):
    s_simpson += 4*f(a+k*h)

for k in range(2, N, 2):
    s_simpson += 2*f(a+k*h)

print("Part A.i - Simpson's rule : ", np.exp((-1)*(b**2))* (1/3) * h * s_simpson) # Equation 5.9

# Scipy.special.dawsn

print("Part A.i - Scipy.special.dawsn : ", special.dawsn(b))
print("====================================================")

"""
Q2 Part A.ii
Purpose - Determine the number of slices needed for an error of O(10^-19) for each method

# Pseudocode
# Repeat Part A.i, but with a different number of slices
# Subtract the approximation from the exact value
# Adjust the number of slices until the error is O(10^-19)
"""

N = 100000                              # Number of slices
a = 0.0                                 # Start Point
b = 4.0                                 # End Point
h = (b-a)/N                             # Width of a slice

# Trapezoidal rule

start = time.time()                     # Measure start time

s_trapezoidal = 0.5* f(a) + 0.5 * f(b)
for k in range(1, N):
    s_trapezoidal += f(a+k*h)

end = time.time()                       # Measure end time

print("Part A.ii - Error of Trapezoidal rule : ", np.exp((-1)*(b**2)) * h * s_trapezoidal - special.dawsn(b))
print("Part A.ii - Number of slices - Trapezoidal rule : ", N)
print("Part A.ii - Timing of Trapezoidal rule's error : ", end - start, "s")

# Simpson's Rule

N = 930                                 # Number of slices
h = (b-a)/N                             # Width of a slice

start = time.time()                     # Measure start time

s_simpson = f(a) + f(b)

for k in range(1, N, 2):
    s_simpson += 4*f(a+k*h)

for k in range(2, N, 2):
    s_simpson += 2*f(a+k*h)

end = time.time()                       # Measure end time

print("Part A.ii - Error of Simpson's rule : ", np.exp((-1)*(b**2))* (1/3) * h * s_simpson - special.dawsn(b))
print("Part A.ii - Number of slices - Simpson's rule : ", N)
print("Part A.ii - Timing of Simpson's rule error : ", end - start, "s")
print("====================================================")

"""
Q2 Part A.iii
Purpose - Demonstrate how to get the error on the estimate without knowing the exact value

# Pseudocode
# Define the number of slices for the first estimate to 32
# Define the number of slices for the second estimate to 64
# Use equation 5.28 and 5.29 from the textbook to get the error on the second estimate
"""

# Trapezoidal rule
N1_trapezoidal = 32                             # Number of slices for the first estimate
N2_trapezoidal = 2 * N1_trapezoidal             # Number of slices for the second estimate
h1_trapezoidal = (b-a)/N1_trapezoidal           # Width of a slice for the first estimate
h2_trapezoidal = (b-a)/N2_trapezoidal           # Width of a slice for the second estimate

s1_trapezoidal = 0.5 * f(a) + 0.5 * f(b)
s2_trapezoidal = 0.5 * f(a) + 0.5 * f(b)

for k in range(1, N1_trapezoidal):
    s1_trapezoidal += f(a+k*h1_trapezoidal)

for k in range(1, N2_trapezoidal):
    s2_trapezoidal += f(a+k*h2_trapezoidal)

I1_trapezoidal = np.exp((-1)*(b**2)) * h1_trapezoidal * s1_trapezoidal
I2_trapezoidal = np.exp((-1)*(b**2)) * h2_trapezoidal * s2_trapezoidal

print("Part A.iii - Error of Second estimate - Trapezoidal rule: ", (1/3) * (I2_trapezoidal - I1_trapezoidal)) # Equation 5.28


# Simpson's Rule
N1_simpson = 32                                 # Number of slices for the first estimate
N2_simpson = 2 * N1_simpson                     # Number of slices for the second estimate
h1_simpson = (b-a)/N1_simpson                   # Width of a slice for the first estimate
h2_simpson = (b-a)/N2_simpson                   # Width of a slice for the second estimate

s1_simpson = f(a) + f(b)
s2_simpson = f(a) + f(b)

for k in range(1, N1_simpson, 2):
    s1_simpson += 4*f(a+k*h1_simpson)

for k in range(1, N2_simpson, 2):
    s2_simpson += 4*f(a+k*h2_simpson)

for k in range(2, N1_simpson, 2):
    s1_simpson += 2*f(a+k*h1_simpson)

for k in range(2, N1_simpson, 2):
    s2_simpson += 2*f(a+k*h2_simpson)

I1_simpson = np.exp((-1)*(b**2))* (1/3) * h1_simpson * s1_simpson
I2_simpson = np.exp((-1)*(b**2))* (1/3) * h2_simpson * s2_simpson

print("Part A.iii - Error of Second estimate - Simpson's rule : ", (1/15) * (I2_simpson - I1_simpson)) # Equation 5.29
print("====================================================")

"""
Q2 Part B
Purpose - Reproduce the diffraction pattern of light

# Pseudocode
# Define a function for the bessel function based on Simpson's rule
# Create a plot of the bessel function for three different m values (m=0,1,2)
# Define a function for the intensity of the light
# Define the wavelength and calculate k
# Plot the density plot of the intensity of the light
"""

# Exercise 5.4 a

def J(m, x):                                        # Bessel function based on Simpson's rule
    N = 1000
    a = 0
    b= np.pi
    h = (b-a)/N

    def g(theta):
        return np.cos(m*theta - x*np.sin(theta))

    s_j = g(a) + g(b)

    for k in range(1, N, 2):
        s_j += 4 * g(a+k*h)
    for k in range(2, N, 2):
        s_j += 2 * g(a+k*h)

    return (1/np.pi) * (1/3) * h * s_j

x_array = np.arange(0, 20.1, 0.1)                     # Create an array from x=0 to x=20

plt.figure(1)                                         # Plot
plt.plot(x_array, J(0, x_array), label="m = 0")
plt.plot(x_array, J(1, x_array), label="m = 1")
plt.plot(x_array, J(2, x_array), label="m = 2")

plt.xlabel("x")
plt.ylabel("J(x)")
plt.title("Bessel Function - Simpson's Rule (m = 0,1,2 and x-interval = 0.1)")
plt.legend()


plt.figure(2)
plt.plot(x_array, special.jv(0, x_array), label="m = 0")
plt.plot(x_array, special.jv(1, x_array), label="m = 1")
plt.plot(x_array, special.jv(2, x_array), label="m = 2")
plt.xlabel("x")
plt.ylabel("J(x)")
plt.title("Bessel Function - Scipy.special.jv (m=0,1,2 and x-interval=0.1)")
plt.legend()




# Exercise 5.4.b


delta = 0.5*(10**(-8))                              # interval of the focal plane
x = y = np.arange(-10**(-6), 10**(-6), delta)
X, Y = np.meshgrid(x, y)

wavelength = 500 * (10**(-9)) # Wavelength : 500 nm     # Wavelength
k = 2*np.pi/wavelength # k = 2pi/lambda


def intensity(kr):                                  # Intensity of the light
    return (J(1, kr) / kr)**2



Z = np.array(intensity(k*np.sqrt(X**2+Y**2)))



plt.figure(3)
plt.imshow(Z, cmap="hot", vmax=0.001, extent=[-1, 1, -1, 1])
plt.xlabel("x (10^-6 meter)")
plt.ylabel("y (10^-6 meter)")
plt.title("Density plot of intensity of the light (Wavelength = 500nm)")
plt.show()


