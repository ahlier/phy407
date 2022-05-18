import numpy as np
import matplotlib.pyplot as plt
import random as rd


''' 
This program will return the square root of mean square displacement vs time

The code is similar to question 2, but without boundary condition, and for each 
step, the displacement of particle is calculated and recorded.
'''

XMAX = YMAX = ZMAX = 10.

dx = dy = dz = 0.1

D = 1

nx, ny, nz = int(2 * XMAX / dx), int(2 * YMAX / dy), int(2 * ZMAX / dz)

dx2, dy2, dz2 = dx * dx, dy * dy, dz * dz
dt = 1e-5  # seconds
tmax = 0.001  # seconds

u0 = np.zeros((nx, ny, nz), int)

u0[int(nx / 2), int(ny / 2), int(nz / 2)] = 1000  # u0[100, 100, 100] = 100
u_array = [u0.copy()]
t_array = [0]


# Function to generate a random number between 0, 1, 2, 3, 4, 5
def random_generator():  # 0 - x up, 1 - x down, 2 - y up, 3 - y down, 4 - z up, 5 - z down
    return rd.randrange(7)


displacement = []  # used to stored displacement for every particle
mean_square_d = []  # used to stored average displacement for particle
t = 0
while t < tmax:
    displacement = []
    u = np.zeros((nx, ny, nz), int)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for i in range(u0[x, y, z]):
                    displacement.append(np.sqrt(
                        (x - 100) ** 2 + (y - 100) ** 2 + (z - 100) ** 2))
                    move = random_generator()
                    new_x = x
                    new_y = y
                    new_z = z

                    if move == 0:
                        new_x += 1
                    elif move == 1:
                        new_x -= 1
                    elif move == 2:
                        new_y += 1
                    elif move == 3:
                        new_y -= 1
                    elif move == 4:
                        new_z += 1
                    elif move == 5:
                        new_z -= 1


                    u[new_x, new_y, new_z] += 1

    mean_square_d.append(np.mean(displacement))

    t += dt
    u0 = u.copy()
    u_array.append(u.copy())
    t_array.append(t)



plt.plot(t_array[1:], mean_square_d)
plt.xlabel('Time')
plt.ylabel('sqrt of MSD')
plt.title('Change in Square Root of Mean Square Displacement \n'
          'with Respect to Change in Time')
plt.show()
