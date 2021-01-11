"""
Purpose - Plot the relief map for Lake Geneva

# Pseudo code
# File Read
#   Find the location of Lake Geneva - 46'44N 6.53E
#   Download the dataset from http://dds.cr.usgs.gov/srtm/version2_1/SRTM3
#   Place the downloaded file into the same folder
#   Import struct, numpy, matplotlib
#   Create a 1201x1201 array to store the height data
#   Read signed 2 byte integers from the file using struct.unpack('>h', buf)[0] and store it to the array
# Calculate Gradient
#   Create a 1201x1201 array to store the gradient data (one for dw/dx, one for dw/dy)
#   Use central difference for interior points, forward difference for the first point, and backward difference for
#   the last point to calculate the derivative
#   Store the calculated value to the 1201x1201 array
# Calculate the illumination
#   Create a 1201x1201 array to store the illumination data
#   Define the angle (pi/6)
#   Define the distance between grid points
#   Use for loop to compute the illumination at each point and store it to the 1201x1201 array
# Plot
#   Use imshow - cmap="gray", extent=(6, 7, 46, 47)
#   Adjust vmin and vmax to get informative plot
#   Plot w and I
"""

import struct
import numpy as np
import matplotlib.pyplot as plt


filename = "N46E006.hgt"
f = open(filename, 'rb')


w = np.zeros(shape=(1201, 1201)) # Create multi-dimensional array 1201 x 1201
dwdx = np.zeros_like(w)
dwdy = np.zeros_like(w)
dx = 83. # distance in meters between grid points
dy = 83. # distance in meters between grid points

for i in range(len(w)):
    for j in range(len(w)):
        buf = f.read(2)
        val = struct.unpack('>h', buf)[0]
        w[i][j] = val

"""
for i in range(len(w)):
    for j in range(len(w)):
        if w[i][j] <= 100.:
            if j == 0:
                w[i][j] = w[i][j+1]
            else:
                w[i][j] = w[i][j-1]
"""


# Compute partial derivative in terms of x
for i in range(len(dwdx)):
    for j in range(len(dwdx)):
        if j == 0:
            dwdx[i][j] = (w[i][j+1] - w[i][j])/dx # Forward difference
        elif j == (len(dwdx)-1):
            dwdx[i][j] = (w[i][j] - w[i][j-1])/dx # Backward difference
        else:
            dwdx[i][j] = (w[i][j+1] - w[i][j-1])/(2*dx) # Central difference

# Compute partial derivative in terms of y
for j in range (len(dwdy)):
    for i in range(len(dwdy)):
        if i == 0:
            dwdy[i][j] = (w[i+1][j] - w[i][j])/dy # Forward difference
        elif i == (len(dwdy)-1):
            dwdy[i][j] = (w[i][j] - w[i-1][j])/dy # Backward difference
        else:
            dwdy[i][j] = (w[i+1][j] - w[i-1][j])/(2*dy) # Central difference



# Compute Illumination
illumination = np.zeros_like(w)
angle = np.pi/6

def compute_illumination(dwdx, dwdy, angle):
    return (-1)*(np.cos(angle)*dwdx + np.sin(angle)*dwdy)/(np.sqrt(dwdx**2 + dwdy**2 + 1))

for i in range(len(illumination)):
    for j in range(len(illumination)):
        illumination[i][j] = compute_illumination(dwdx[i][j], dwdy[i][j], angle)





# Plot

plt.figure(1) # largest = 3145
plt.imshow(w, cmap="gist_earth", vmin=0, vmax = 3000, extent=(6, 7, 46, 47))
plt.xlabel("Longitude (in degree)")
plt.ylabel("Latitude (in degree)")
plt.title("Height of Lake Geneva - $46\degree$26'N $6\degree$33'E")
plt.colorbar()

plt.figure(2) # largest = 422
plt.imshow(dwdx, cmap="gray", vmin=-0.5, vmax=0.5, extent=(6, 7, 46, 47))
plt.xlabel("Longitude (in degree)")
plt.ylabel("Latitude (in degree)")
plt.title("dwdx of Lake Geneva - $46\degree$26'N $6\degree$33'E")
plt.colorbar()

plt.figure(3)
plt.imshow(dwdy, cmap="gray", vmin=-0.5, vmax=0.5, extent=(6, 7, 46, 47))
plt.xlabel("Longitude (in degree)")
plt.ylabel("Latitude (in degree)")
plt.title("dwdy of Lake Geneva - $46\degree$26'N $6\degree$33'E")
plt.colorbar()

plt.figure(4) # largest = -0.403549, smallest = -0.500137
plt.imshow(illumination, cmap="gray", vmin = -0.3, vmax = 0.3, extent=(6, 7, 46, 47))
plt.xlabel("Longitude (in degree)", fontsize=16)
plt.ylabel("Latitude (in degree)", fontsize=16)
plt.title("Illuminated of Lake Geneva - $46\degree$26'N $6\degree$33'E", fontsize=16)
plt.colorbar()


plt.figure(5) # Zoomed: CERN
plt.imshow(w[:601][:601], cmap="gist_earth", vmin=0, vmax = 3000, extent=(6, 6.5, 46, 46.5))
plt.xlabel("Longitude (in degree)")
plt.ylabel("Latitude (in degree)")
plt.title("CERN - $46\degree$23'N $6\degree$05'E")
plt.colorbar()

plt.figure(6) # largest = -0.403549, smallest = -0.500137
plt.imshow(illumination[:601][:601], cmap="gray", vmin = -0.3, vmax = 0.3, extent=(6, 6.5, 46, 46.5))
plt.xlabel("Longitude (in degree)", fontsize=16)
plt.ylabel("Latitude (in degree)", fontsize=16)
plt.title("CERN - $46\degree$23'N $6\degree$05'E", fontsize=16)
plt.colorbar()

plt.show()
f.close()