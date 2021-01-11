import numpy as np
import matplotlib.pyplot as plt
import time
"""
Question 1 - Part a
Purpose - Create a contour plot (density plot) of the potential and a stream plot of the electric field lines, 
          using Gauss-Seidel Method without overrelaxation

# Pseudocode

# Potential
# Create constnats
#   M: Number of grid points on each side
#   V: 1.0 volt
#   target: Target accuracy
#   a: spacing between grid points in meter
#   ftsz: fontsize
# Create an array to hold the values of potential
#   phi: 100 x 100 array
# Set the voltage to 1V and -1V to two plates as shown in Figure 3
# Perform a while loop
#   Condition: error > target accuracy 
#   Perform a nested for-loop
#       Excluding two plates and the boundaries, set the potential at i,j based on Equation 9.16 from the textbook
#       Compute the difference between the old and new value of the potential and check against the target accuracy
# Plot a contour plot

# Electric Field
# Create arrays to hold the values of electric field
#   Ex: Electric field in x direction
#   Ey: Electric field in y direction
# Compute the electric field by taking a derivative of the potential (i.e E = -V')
#   Forward difference: for the left-most position
#   Backward difference: for the right-most position
#   Central difference: for inner position
# Plot a stream plot with the values of electric field
"""
# Constants
M = 100             # Grid Squares
V = 1.0             # Voltage on a plate
target =1e-6        # Target accuracy
a = 0.001            # 0.1 cm = 0.001 m, spacing between grid points
ftsz = 20


# Create an array to hold potential values
phi = np.zeros([M, M], float)           # 101 x 101 array
phi[20:81, 20] = V                      # First plate with 1V voltage
phi[20:81, 80] = -1 * V                 # Second plate with -1V voltage

start_time1 = time.time()

# Main Loop
delta = 2 * target
while delta > target:
    delta = 0.0
    for i in range(M):
        for j in range(M):
            if not (i == 0) and not (i == M-1) and not (j == 0) and not (j == M-1) and not (20 <= i <= 80 and (j == 20 or j == 80)):
                phi_old = phi[i, j]
                phi_new = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1]) / 4
                phi[i, j] = phi_new

                difference = abs(phi_new - phi_old) # Find the difference between old and new potential

                if difference > delta: # Overwrite delta if the difference is bigger than the value of deltat that's currently stored
                    delta = difference
end_time1 = time.time()

# Density Plot (No Overrelaxation)
plt.figure(1)
plt.imshow(phi, extent=[0, 10, 0, 10])
cbar = plt.colorbar()
cbar.set_label("Potential (in V)", fontsize=ftsz)
cbar.ax.tick_params(labelsize=16)
plt.title("The Density plot of the electric potential without overrelaxation", fontsize=ftsz)
plt.xlabel("$x$ (in cm)", fontsize=ftsz)
plt.ylabel("$y$ (in cm)", fontsize=ftsz)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)


# Calculation of Electric field
Ex = np.zeros([M, M], float)    # Electric field in x direction
Ey = np.zeros([M, M], float)    # Electric field in y direction
x = np.linspace(0, 10, 100)
y = np.linspace(0, 10, 100)
X, Y = np.meshgrid(x, y)

for i in range(M):
    for j in range(M):
        if j == 0: # Forward Difference
            Ex[i, j] = (-1) * (phi[i, j+1] - phi[i, j]) / a
        elif j == M-1: # Backward Difference
            Ex[i, j] = (-1) * (phi[i, j] - phi[i, j-1]) / a
        else:       # Central Difference
            Ex[i, j] = (-1) * (phi[i, j+1] - phi[i, j-1]) / a


        if i == 0: # Forward Difference
            Ey[i, j] = (-1) * (phi[i+1, j] - phi[i, j]) / a
        elif i == M-1: # Backward Difference
            Ey[i, j] = (-1) * (phi[i, j] - phi[i-1, j]) / a
        else:           # Central Difference
            Ey[i, j] = (-1) * (phi[i+1, j] - phi[i-1, j]) / a

# Plot stream-plot
fig = plt.figure(2, figsize=(6, 3))
strm = plt.streamplot(X, Y, Ex, Ey, color=phi, linewidth=2, cmap='afmhot') #afmhot
cbar = fig.colorbar(strm.lines)
cbar.set_label('Potential (in V)', fontsize=ftsz)
cbar.ax.tick_params(labelsize=16)
plt.title('Electric Field Lines', fontsize=ftsz)
plt.xlabel("$x$ (in cm)", fontsize=ftsz)
plt.ylabel("$y$ (in cm)", fontsize=ftsz)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

# Plot contour-plot
fig = plt.figure(3, figsize=(8, 6))
ctf = plt.contourf(X, Y, phi, 50, cmap="afmhot")
cbar = plt.colorbar()
cbar.set_label("Potential (in V)", fontsize=ftsz)
cbar.ax.tick_params(labelsize=16)
plt.title("The Contour plot of the electric potential without overrelaxation", fontsize=ftsz)
plt.xlabel("$x$ (in cm)", fontsize=ftsz)
plt.ylabel("$y$ (in cm)", fontsize=ftsz)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)


"""
Question 1 - Part b
Purpose - Repeat Part a with Gauss-Seidel method with overrelaxation

# Pseudocode
# Set Constants
#   omega: overrelaxation parameter
# Create an array to hold the values of potential
#   phi: 100 x 100 array
# Set the voltage to 1V and -1V to two plates as shown in Figure 3
# Perform a while loop
#   Condition: error > target accuracy 
#   Perform a nested for-loop
#       Excluding two plates and the boundaries, set the potential at i,j based on Equation 9.17 from the textbook
#       Compute the difference between the old and new value of the potential and check against the target accuracy
# Plot a contour plot
"""

# Constant
omega = 0.1

# Rewrite the array
phi = np.zeros([M, M], float)           # 101 x 101 array
phi[20:81, 20] = V                      # First plate with 1V voltage
phi[20:81, 80] = -1 * V                 # Second plate with -1V voltage

start_time2 = time.time()
# Main Loop
delta = 2 * target
while delta > target:
    delta = 0.0
    for i in range(M):
        for j in range(M):
            if not (i == 0) and not (i == M-1) and not (j == 0) and not (j == M-1) and not (20 <= i <= 80 and (j == 20 or j == 80)):
                phi_old = phi[i, j]
                phi_new = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1]) * (1 + omega) / 4 - (omega * phi_old)
                phi[i, j] = phi_new

                difference = abs(phi_new - phi_old) # Find the difference between old and new potential

                if difference > delta: # Overwrite delta if the difference is bigger than the value of deltat that's currently stored
                    delta = difference
end_time2 = time.time()

# Contour Plot (Overrelaxation - omage = 0.1)
plt.figure(4)
plt.imshow(phi, extent=[0, 10, 0, 10])
cbar = plt.colorbar()
cbar.set_label("Potential (in V)", fontsize=ftsz)
cbar.ax.tick_params(labelsize=16)
plt.title("The contour plot of the electric potential with over-relaxation (w = 0.1)", fontsize=ftsz)
plt.xlabel("$x$ (in cm)", fontsize=ftsz)
plt.ylabel("$y$ (in cm)", fontsize=ftsz)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)


# Constant
omega = 0.5

# Rewrite the array
phi = np.zeros([M, M], float)       # 100 x 100 array
phi[20:81, 20] = V                      # First plate with 1V voltage
phi[20:81, 80] = -1 * V                 # Second plate with -1V voltage

start_time3 = time.time()

# Main Loop
delta = 2 * target
while delta > target:
    delta = 0.0
    for i in range(M):
        for j in range(M):
            if not (i == 0) and not (i == M-1) and not (j == 0) and not (j == M-1) and not (20 <= i <= 80 and (j == 20 or j == 80)):
                phi_old = phi[i, j]
                phi_new = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1]) * (1 + omega) / 4 - (omega * phi_old)
                phi[i, j] = phi_new

                difference = abs(phi_new - phi_old) # Find the difference between old and new potential

                if difference > delta: # Overwrite delta if the difference is bigger than the value of deltat that's currently stored
                    delta = difference
end_time3 = time.time()

# Contour Plot (Overrelaxation - omage = 0.9)
plt.figure(5)
plt.imshow(phi, extent=[0, 10, 0, 10])
cbar = plt.colorbar()
cbar.set_label("Potential (in V)", fontsize=ftsz)
cbar.ax.tick_params(labelsize=16)
plt.title("The contour plot of the electric potential with over-relaxation (w = 0.9)", fontsize=ftsz)
plt.xlabel("$x$ (in cm)", fontsize=ftsz)
plt.ylabel("$y$ (in cm)", fontsize=ftsz)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)


# Print the run time of each case

print("The execution time of Gauss-Seidel method without overrelaxation is", end_time1 - start_time1,"s")
print("The execution time of Gauss-Seidel method with overrelaxation (w = 0.1) is", end_time2 - start_time2, "s")
print("The execution time of Gauss-Seidel method with overrelaxation (w = 0.5) is", end_time3 - start_time3, "s")
plt.show()