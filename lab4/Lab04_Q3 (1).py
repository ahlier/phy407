import numpy as np
import matplotlib.pyplot as plt
from math import exp

"""
Exercise 6.10 Part a

# Pseudocode
# Set the initial guess
# Set the known parameter c
# Perform for loop to implement relaxation method

x=1 converges to 0.7968121300200202

"""

c = 2
x = 1.

print("Exercise 6.10 Part a")
for k in range(50):
    x = 1 - exp((-1)*c*x)
    print("Exercise 6.10 Part a - ", k+1, "th iteration:",x)

print("Note that the value of x converges to", x)
print("======================================================================================================")


"""
Exercise 6.10 Part b
Purpose - Plot c vs x

# Pseudocode
# Create a numpy array for 'c'
# Create a numpy array for 'x' with the length as the array 'c'
# Perform for loop to store the converged value of x
# Plot using pyplot
"""

c = np.linspace(0, 3, 300)
x = np.ones(len(c))

for i in range(len(c)):
    for k in range(50):
        x[i] = 1 - exp((-1)*c[i]*x[i])

plt.figure(1)
plt.errorbar(c, x, color='red', fmt='.', markersize=3, label="x = 1 - exp(-cx)")
plt.xlabel("c (Known parameter)", fontsize=20)
plt.ylabel("x (Unknown parameter)", fontsize=20)
plt.title("c vs x from x = 1 - exp(-cx) by relaxation method", fontsize=20)
plt.legend(fontsize=20)


"""
Exercise 6.11 Part b
Purpose - Print out the number of iterations taken to converge to a solution accurate to 10^-6, using relaxation method

# Set the parameter 'c'
# Set a counter for iteration
# Define a function for the error of relaxation method (Eq6.87)
# Perform for loop to implement eq 6.87
"""


c = 2
iteration = 0
x0 = 1. # Starting point stored for output
x = 1. # Starting point

def error_relaxation(x, x_prime, x_prime2):
    return (x_prime - x_prime2)**2/(2*x_prime - x - x_prime2)

for k in range(50):
    x_prime = 1 - exp((-1)*c*x)
    x_prime2 = 1 - exp((-1)*c*x_prime)


    iteration += 1
    error = error_relaxation(x, x_prime, x_prime2)
    if abs(error) < 10**(-6):
        print("Exercise 6.11 Part b - Relaxation Method")
        print("Exercise 6.11 Part b - The estimate (starting point x =", x0, ", c =", c, ") is", x_prime2)
        print("Exercise 6.11 Part b - Starting from x =", x0,", the error on", iteration+1, "th iteration is", error)
        print("This is close to a solution with an accuracy of less than 10^-6")
        break

    x = x_prime

print("======================================================================================================")

"""
Exercise 6.11 Part c
Purpose - Repeat part b, but using overrelaxation method

# Pseudocode
# Create an array of relaxation factor 'w'
# Define a function for the error of overrelaxation
# Perform for loop to implement the equation in ex6.11
"""


c = 2
iteration = 0
x0 = 1.
x = 1.
w = np.array([0.5, 1, 1.5, 2])  # list of relaxation parameter

def error_overrelaxation(x, x_prime):
    return (x - x_prime)/(1 - 1/((1+w)*(c*exp(-1*c*x)) - w))

for w in w:
    for k in range(50):

        x_prime = (1 + w)*(1 - exp((-1) * c * x)) - w * x

        iteration += 1
        error = error_overrelaxation(x, x_prime)
        if abs(error) < 10**(-6):
            print("Exercise 6.11 Part c - Overrelaxation Method")
            print("Exercise 6.11 Part c - The estimate (w =", w, ", starting point x =", x0, ", c =", c, ") is", x_prime)
            print("Exercise 6.11 Part c - Starting from x =", x0, ", the error on", iteration, "th iteration is", error)
            print("This is close to a solution with an accuracy of less than 10^-6")
            print("======================================================================================================")
            break

        x = x_prime

"""
Exercise 6.13 Part b
Purpose - Solve the Wien displacement law using three methods: Binary search, Newton's, relaxation

# Pseudocode
# Set Planck's constant, speed of light, and Boltzmann constant
# Binary Search
#   Set initial two points for an interval
#   Set target accuracy
#   Set a counter for iteration
#   Perform while loop to implement the four steps in the textbook for binary search
# Relaxation Method 
#   Repeat Ex 6.11 part b
# Newton's Method
#   Set initial guess
#   Set a counter for iteration
#   Use for loop to implement equation 6.96 

x converges to 4.965114231744276 
"""

# Binary Search
h = 6.62607004 * 10**(-34)          # Planck's constant, Unit: m^2 * kg / s
c = 3* 10**8                        # Speed of Light, Unit: m/s
k_b = 1.38064852 * 10**(-23)        # Boltzmann's constant, Unit:m^2 * kg * s^-2 * K^-1

def wien(x):
    return 5*exp((-1)*x) + x - 5

x1 = 2
x2 = 8
target_accuracy = 10**(-6)
iteration = 0

while abs(x1 - x2) > target_accuracy:
    if wien(x1) * wien(x2) < 0:
        mid_point = (x1 + x2)/2
        if wien(mid_point) * wien(x1) > 0:
            x1 = mid_point
        else:
            x2 = mid_point
    iteration += 1

wien_constant_binary = h*c/(k_b*0.5*(x1 + x2))

print("Exercise 6.13 Part b - Binary Search")
print("Exercise 6.13 Part b - The final estimate of the position of the root is", 0.5*(x1 + x2))
print("Exercise 6.13 Part b - The relative error is", abs(x1 - x2))
print("Exercise 6.13 Part b - The number of iterations is", iteration)
print("Exercise 6.13 Part b - The Wien displacement constant is", wien_constant_binary)

print("======================================================================================================")

# Relaxation Method
x0 = 8.
x = 8.
iteration = 0

def error_relaxation(x, x_prime, x_prime2):
    return (x_prime - x_prime2)**2/(2*x_prime - x - x_prime2)

for k in range(50):
    x_prime = 5 - 5*exp((-1)*x)
    x_prime2 = 5 - 5*exp((-1)*x_prime)


    iteration += 1
    error = error_relaxation(x, x_prime, x_prime2)
    if abs(error) < 10**(-6):
        print("Exercise 6.13 Part b - Relaxation Method")
        print("Exercise 6.13 Part b - The position (starting point x =", x0, ") is", x_prime2)
        print("Exercise 6.13 Part b - Starting from x =", x0,", the error on", iteration+1, "th iteration is", error)
        print("This is close to a solution with an accuracy of less than 10^-6")
        print("The Wien displacement constant is", h*c/(k_b*x_prime2))
        wien_constant_relaxation = h*c/(k_b*x_prime2)
        break

    x = x_prime


print("======================================================================================================")

# Newton's Method
x0 = 8.
x = 8.
iteration = 0

for k in range(50):
    x_prime = x - (5*exp(-1*x) + x - 5)/(-5*exp(-1*x)+1)

    iteration += 1
    error = x_prime - x

    if abs(error) < 10**(-6):
        print("Exercise 6.13 Part b - Newton's method")
        print("Exercise 6.13 Part b - The position (starting point x =", x0, ") is", x_prime)
        print("Exercise 6.13 Part b - Starting from x =", x0,", the error on", iteration, "th iteration is", error)
        print("This is close to a solution with an accuracy of less than 10^-6")
        print("The Wien displacement constant is", h*c/(k_b*x_prime))
        wien_constant_newton = h*c/(k_b*x_prime)
        break

    x = x_prime

print("======================================================================================================")

"""
Exercise 6.13 Part c
Purpose - Calculate the surface temperature of the Sun with Wien constant

Pseudocode
# Divide Wien constant by the given wavelength
"""

wavelength = 502* 10**(-9) # Unit: m
print("Exercise 6.13 Part c - Binary Search - The surface temperature of the sun is", wien_constant_binary/wavelength)
print("Exercise 6.13 Part c - Relaxation Method - The surface temperature of the sun is", wien_constant_relaxation/wavelength)
print("Exercise 6.13 Part c - Newton's Method - The surface temperature of the sun is", wien_constant_newton/wavelength)
plt.show()