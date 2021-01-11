import numpy as np
import matplotlib.pyplot as plt


"""
Q2 - Part a
Purpose : Create a pseudocode of a function that computes the equation (12)

# Pseudocode
# As code comments
# From keyboard, read the maximum population, r
# From keyboard, read the number of years, p
# From keyboard, read the initial value of x, x0
# Create an array with length of p
# Compute x_p using equation 12 and assign it to the array
# Plot the array
"""


"""
Q2 - Part b
Purpose : Create a function based on the pseudocode from part a
"""
def compute_array(x0, r, p_max):
    x_array = np.zeros(p_max)
    x_array[0] = x0

    x_temp = x0
    for i in range(1, p_max):
        x_temp = r*(1-x_temp)*x_temp
        x_array[i] = x_temp

    return x_array

"""
Q2 - Part c
Purpose : Create a plot of population growth at different maximum reproduction rate from 2.0 to 4.0 

# Pseudocode
# Define constants - initial population & the number of years
# Create an array for years
# Create an array for population by using the function defined in part b
# Plot
"""


x0 = 0.1                                                    # Initial population = 0.1
p_max = 51                                                  # Number of years = 50

p_array = np.arange(p_max)
x_array1 = compute_array(x0, 2, p_max)                      # The second parameter represents the maximum reproduction rate
x_array2 = compute_array(x0, 2.2, p_max)
x_array3 = compute_array(x0, 2.4, p_max)
x_array4 = compute_array(x0, 2.6, p_max)
x_array5 = compute_array(x0, 2.8, p_max)
x_array6 = compute_array(x0, 3.0, p_max)
x_array7 = compute_array(x0, 3.2, p_max)
x_array8 = compute_array(x0, 3.4, p_max)
x_array9 = compute_array(x0, 3.6, p_max)
x_array10 = compute_array(x0, 3.8, p_max)
x_array11 = compute_array(x0, 4.0, p_max)

plt.figure(1, figsize=(8,6))                                # Plot
plt.plot(p_array, x_array1, label = "x0 = 0.1, r = 2")
plt.plot(p_array, x_array2, label = "x0 = 0.1, r = 2.2")
plt.plot(p_array, x_array3, label = "x0 = 0.1, r = 2.4")
plt.plot(p_array, x_array4, label = "x0 = 0.1, r = 2.6")
plt.plot(p_array, x_array5, label = "x0 = 0.1, r = 2.8")
plt.plot(p_array, x_array6, label = "x0 = 0.1, r = 3.0")
#plt.plot(p_array, x_array7, label = "x0 = 0.1, r = 3.2")
#plt.plot(p_array, x_array8, label = "x0 = 0.1, r = 3.4")
plt.plot(p_array, x_array9, label = "x0 = 0.1, r = 3.6")
#plt.plot(p_array, x_array10, label = "x0 = 0.1, r = 3.8")
plt.plot(p_array, x_array11, label = "x0 = 0.1, r = 4.0")
plt.xlabel("Years")
plt.ylabel("Population")
plt.title("Years vs Population")
plt.legend()



"""
Q2 - Part d
Purpose : Create a plot of population growth at different maximum reproduction rate from 2.0 to 4.0 with the increment of 0.1
Double Period - After 3
chaos - After 3.570

# Pseudocode
# Create an array for maximum reproduction rates - from 2 to 4, increment: 0.1
# Define the number of years
# Create an array for years
# Create an array for population at each rate
# Plot
"""

r_array = np.arange(2, 4.1, 0.1)                             # Increment of maximum reproduction rate : 0.1
for i in range(len(r_array)):                                # Cut the rate to the first decimal (Float)
    r_array[i] = round(r_array[i], 1)

p_max = 2001                                                 # Number of years = 2000
p_array = np.arange(p_max)

part_d_dict = {}                                             # Store the population array at each rate in a dictionary
for r in r_array:
    part_d_dict[r] = compute_array(x0, r, p_max)


plt.figure(2, figsize=(8, 6))                                # Plot
for r in r_array:
    if r < 3:
        plt.plot([r] * len(part_d_dict[r][-100:]), part_d_dict[r][-100:], 'k.', markersize=0.1)
    else:
        plt.plot([r] * len(part_d_dict[r][-1000:]), part_d_dict[r][-1000:], 'k.', markersize=0.1)
plt.xlabel("Maximum Reproduction Rate")
plt.ylabel("Population")
plt.title("Maximum Reproduction Rate vs Population - Increment of the rate : 0.1")


r_array2 = np.arange(2, 4.005, 0.005)                       # Increment of maximum reproduction rate : 0.005
for i in range(len(r_array2)):                              # Cut the rate to the third decimal (Float)
    r_array2[i] = round(r_array2[i], 3)

part_d_dict2 = {}                                           # Store the population array at each rate in a dictionary
for r in r_array2:
    part_d_dict2[r] = compute_array(x0, r, p_max)

plt.figure(3, figsize=(8, 6))                               # Plot - Last 100 items for r < 3 and last 1000 items otherwise
for r in r_array2:
    if r < 3:
        plt.plot([r] * len(part_d_dict2[r][-100:]), part_d_dict2[r][-100:], 'k.', markersize=0.1)
    else:
        plt.plot([r] * len(part_d_dict2[r][-1000:]), part_d_dict2[r][-1000:], 'k.', markersize=0.1)
plt.xlabel("Maximum Reproduction Rate")
plt.ylabel("Population")
plt.title("Maximum Reproduction Rate vs Population - Increment of the rate : 0.005")

"""
Q2 - Part e
Purpose : Create a plot of population growth at different maximum reproduction rate from 3.738 to 3.745 with the increment of 10^(-5)

# Pseudocode
# Create an array for maximum reproduction rates - from 3.738 to 3.745, increment: 10**(-5)
# Use the same years from part d
# Create an array for population at each rate
# Plot
"""

r_array3 = np.arange(3.738, 3.745, 10**(-5))                # Increment of maximum reproduction rate : 10^(-5)
part_e_dict = {}                                            # Store the population array at each rate in a dictionary
for r in r_array3:
    part_e_dict[r] = compute_array(x0, r, p_max)


plt.figure(4, figsize=(10, 8))                              # Plot - Last 100 items for r < 3.741 and last 1000 items otherwise
for r in r_array3:
    if r < 3.741:
        plt.plot([r] * len(part_e_dict[r][-100:]), part_e_dict[r][-100:], 'k.', markersize=0.1)
    else:
        plt.plot([r] * len(part_e_dict[r][-1000:]), part_e_dict[r][-1000:], 'k.', markersize=0.1)
plt.xlabel("Maximum Reproduction Rate")
plt.ylabel("Population")
plt.title("Maximum Reproduction Rate vs Population - Increment of the rate : 10^(-5)")



"""
Q2 - Part f
Purpose - Set the initial populations and the maximum reproduction rate that show a chaotic behavior

# Pseudocode
# Define the initial population and the difference between two populations
# Define the number of years
# Create an array for years
# Create an array for population at the maximum reproduction rate of 4.0
# Plot
"""

e = x0*(10**(-7))                                           # x0 = 0.1
x1 = x0                                                     # e = difference in the initial population
x2 = x0 + e                                                 # x1, x2 = initial populations
p_max = 201                                                 # Number of years = 200
p_array = np.arange(p_max)
x_array_e1 = compute_array(x1, 4.0, p_max)                  # Create a population array with maximum reproduction rate of 4.0
x_array_e2 = compute_array(x2, 4.0, p_max)

plt.figure(5, figsize=(8, 6))                               # Plot
plt.plot(p_array, x_array_e1, label = "x1 = 0.1, r = 4.0")
plt.plot(p_array, x_array_e2, label = "x2 = 0.100001, r = 4.0")
plt.xlabel("Years")
plt.ylabel("Population")
plt.title("Years vs Population - Q2.f")
plt.legend()



"""
Q2 - Part g
Purpose : Plot the difference between two population histories in a logarithmic scale and obtain Lypunov exponent by curve_fitting
"""
# Pseudocode
# Take the absolute value of differences between two population histories from part f
# Define a model function : Ae^(a*lyapunov)
# Plot - the y-axis in a logarithmic scale
x_array_diff = np.abs(x_array_e2 - x_array_e1)              # Absolute difference of two population histories

def model_func(p, a, lyapunov):                             # Model function : Lyapunov exponent
    return a * np.exp(lyapunov * p)

plt.figure(6, figsize=(8, 6))                               # Plot - y-axis in a logarithmic scale
plt.semilogy(p_array, x_array_diff)                         # Fit using the model function
plt.semilogy(p_array[:26], model_func(p_array[:26], x_array_diff[0], 0.7), label="Best Fit")
plt.xlabel("Years")
plt.ylabel("log(Population Difference)")
plt.title("Years vs log(Population Difference) - Q2.f")
plt.legend()

plt.show()