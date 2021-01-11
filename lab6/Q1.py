import numpy as np
import matplotlib.pyplot as plt

'''
This program simulate the situation where a ball is orbiting around a rod due
to solely on gravitational force

pseudo code:
- write out function of force
- run a for-loop, and for each iteration, apply the Runge-Kutta method to get
the x, y position, and store them
- make a plot using x,y position
'''
G = 1
M = 10
L = 2
x_i, y_i = 1, 0  # initial x,y position
vx_i, vy_i = 0, 1  # initial x,y velocity
duration = 10  # total time
dt = 0.01  # time interval between each data point

# following code used textbook code page 345 as reference.
def f(data):
    # this function will return position abd velocity (in order of x velocity,
    # x acceleration, y velocity, y acceleration) of object with given position
    # and velocity using given equations
    x = data[0]
    vx = data[1]
    y = data[2]
    vy = data[3]
    r = (x ** 2 + y ** 2)**0.5
    return np.array([vx, -G * M * x / (r ** 2 * (r ** 2 + L ** 2 / 4)**0.5),
                  vy, -G * M * y / (r ** 2 * (r ** 2 + L ** 2 / 4)**0.5)], float)

time = np.arange(0, duration, dt)
xpoints = []  # stored all the x position
ypoints = []  # stored all the y position
data = np.array([x_i, vx_i, y_i, vy_i], float)
for t in time:
    # each iteration will apply the Runge-Kutta method, and store the updated
    # x,y position
    xpoints.append(data[0])
    ypoints.append(data[2])
    k1 = dt * f(data)
    k2 = dt * f(data + 0.5 * k1)
    k3 = dt * f(data + 0.5 * k2)
    k4 = dt * f(data + k3)
    data += (k1 + 2 * k2 + 2 * k3 + k4) / 6

plt.plot(xpoints, ypoints, label='ball')
plt.xlabel('x position')
plt.ylabel('y position')
plt.title('Orbit of a ball around a rod')
plt.legend()
plt.show()
