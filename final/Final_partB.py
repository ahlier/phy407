import numpy as np
import matplotlib.pyplot as plt
import random as rd

XMAX = YMAX = ZMAX = 10.

dx = dy = dz = 0.1

D = 1



nx, ny, nz = int(2*XMAX/dx), int(2*YMAX/dy), int(2*ZMAX/dz)

dx2, dy2, dz2 = dx*dx, dy*dy, dz*dz
dt = 1e-5 # seconds
tmax = 0.001 # seconds

u0 = np.zeros((nx, ny, nz), int)

u0[int(nx/2), int(ny/2), int(nz/2)] = 100000 # u0[100, 100, 100] = 100
u_array = [u0.copy()]
t_array = [0]
# Function to generate a random number between 0, 1, 2, 3, 4, 5
def random_generator():             # 0 - x up, 1 - x down, 2 - y up, 3 - y down, 4 - z up, 5 - z down
    return rd.randrange(6)


# Check if the particle hit the boundary
# Return True if it did
# Return False otherwise
def check_boundary(x, y, z):
    if 0 <= x <= 199 and 0 <= y <= 199 and 0 <= z <= 199:
        return False
    else:
        return True




t = 0
while t < tmax:
    u = np.zeros((nx, ny, nz), int)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for i in range(u0[x, y, z]):
                    outside_boundary = True
                    while(outside_boundary):
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

                        outside_boundary = check_boundary(new_x, new_y, new_z)

                    u[new_x, new_y, new_z] += 1

    t += dt
    u0 = u.copy()
    u_array.append(u.copy())
    t_array.append(t)
    print(u_array[-1][100, 101, 100])
    print("Sum:", u0.sum())
    print("===========================================")

