import numpy as np
import matplotlib.pyplot as plt


N = 2
Lx = 4.0
Ly = 4.0
dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)
x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)
xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()

m = 1  # mass of particle
sigma = 1  # given
epsilon = 1  # given
dt = 0.01  # time interval
T = 1000


def acceleration(x, y):
    a = 4 * epsilon * (12 * sigma ** 12 / (x ** 2 + y ** 2) ** (
            13 / 2) - 6 * sigma ** 6 / (x ** 2 + y ** 2) ** (7 / 2)) / m
    a_x = a * x / np.sqrt(x ** 2 + y ** 2)
    a_y = a * y / np.sqrt(x ** 2 + y ** 2)
    return a_x, a_y


def verlet(dt, x, y, v_x, v_y):
    new_x = x + dt * v_x
    new_y = y + dt * v_y

    a_x, a_y = acceleration(new_x, new_y)
    k_x, k_y = dt * a_x, dt * a_y
    v_new_x = v_x + k_x
    v_new_y = v_y + k_y
    return new_x, new_y, v_new_x, v_new_y


def Q2(p1, p2, v_x, v_y):

    x = p1[0]-p2[0]  # magnitude of x separation
    y = p1[1]-p2[1]  # magnitude of y separation

    p1_x, p1_y = [p1[0]], [p2[1]]  # x,y position of first particle
    p2_x, p2_y = [p1[0]], [p2[1]]

    x, y, v_x, v_y = verlet(dt, x, y, v_x, v_y)

    delta_x = x - (p1_x[-1] - p2_x[-1])  # change in x direction
    delta_y = y - (p1_y[-1] - p2_y[-1])  # change in y direction

    p1_x_new = 0.5 * delta_x + p1_x[-1]
    p1_y_new = 0.5 * delta_y + p1_y[-1]

    return p1_x_new, p1_y_new, v_x, v_y

recorded_x, recorded_y = [], []
recorded_v_x, recorded_v_y = [], []

for i in range(len(x_initial)):
    recorded_x.append([x_initial[i]])
    recorded_y.append([y_initial[i]])
    recorded_v_x.append([0])
    recorded_v_y.append([0])



for i in range(len(x_initial)):
    for k in range(len(x_initial)):
        if k != i:
            x = recorded_x[i][-1] - recorded_x[k][-1]  # magnitude of x separation
            y = recorded_y[i][-1] - recorded_y[k][-1]  # magnitude of y separation

            a_x, a_y = acceleration(x, y)  # x,y component of acceleration
            v_x = 0.5 * dt * a_x  # x velocity in first dt/2 sec
            v_y = 0.5 * dt * a_y  # y velocity in first dt/2 sec

            recorded_v_x[i][-1] += v_x
            recorded_v_y[i][-1] += v_y

for j in range(T):
    x_copy, y_copy = np.copy(recorded_x), np.copy(recorded_y)
    v_x_copy, v_y_copy = np.copy(recorded_v_x), np.copy(recorded_v_y)

    for i in range(len(x_initial)):

        for k in range(len(x_initial)):
            if i != k:
                relative_vx = v_x_copy[i][-1]-v_x_copy[k][-1]
                relative_vy = v_y_copy[i][-1]-v_y_copy[k][-1]
                a, b, c, d = Q2([x_copy[i][-1], y_copy[i][-1]],
                                [x_copy[k][-1], y_copy[k][-1]], relative_vx,
                                relative_vy)
                recorded_x[i][-1] += a/N
                recorded_y[i][-1] += b/N
                recorded_v_x[i][-1] += c/N
                recorded_v_y[i][-1] += d/N

        recorded_v_x[i].append(recorded_v_x[i][-1])
        recorded_v_y[i].append(recorded_v_y[i][-1])
        recorded_x[i].append(recorded_x[i][-1])
        recorded_y[i].append(recorded_y[i][-1])
# plt.plot(recorded_x[0], recorded_y[0])
# plt.show()
print(recorded_x)



