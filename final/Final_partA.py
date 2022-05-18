import numpy as np
import matplotlib.pyplot as plt
import time


XMAX = YMAX = ZMAX = 10.

dx = dy = dz = 0.1

D = 1



nx, ny, nz = int(2*XMAX/dx), int(2*YMAX/dy), int(2*ZMAX/dz)

dx2, dy2, dz2 = dx*dx, dy*dy, dz*dz
dt = 1e-5 # seconds
tmax = 0.0001 # seconds

u0 = np.zeros((nx, ny, nz), int)
u0[int(nx/2), int(ny/2), int(nz/2)] = 10000 # u0[100, 100, 100] = 10000
u = u0.copy()       # Make a local copy
u_array = [u0]
t_array = [0]


def timestep(u0, u):
    u[1:-1, 1:-1, 1:-1] = u0[1:-1, 1:-1, 1:-1] + D * dt * (
            (u0[2:, 1:-1, 1:-1] - 2 * u0[1:-1, 1:-1, 1:-1] + u0[:-2, 1:-1, 1:-1])/dx2
            + (u0[1:-1, 2:, 1:-1] -2 * u0[1:-1, 1:-1, 1:-1] + u0[1:-1, :-2, 1:-1])/dy2
            + (u0[1:-1, 1:-1, 2:] - 2 * u0[1:-1, 1:-1, 1:-1] + u0[1:-1, 1:-1, :-2])/dz2)
    u0 = u.copy()

    return u0, u




t = 0
print("Loop start")
start_time = time.time()
while t < tmax:
    u0, u = timestep(u0, u)
    t += dt
    u_array.append(u.copy())
    t_array.append(t)
    print(t, u[100,100,100])

print("Loop end")
end_time = time.time()
print("Total time:", end_time - start_time, "s")


for i in range(len(t_array)):
    print("Time : ", t_array[i], " seconds with ", u_array[i][100][100][100]," particles, total particles:", u_array[i].sum())

