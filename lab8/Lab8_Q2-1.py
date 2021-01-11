import numpy as np
import matplotlib.pyplot as plt


'''
This program uses 1D FTCS method to simulate the motion of wave.

pseudo code:
-First set the values constants and variables
-Then run a loop, each iteration will use FTCS method update the position and 
speed of the wave, also update the boundary points using forward/backward 
difference method.
-plot the graph in each iteration to create an animation
-make the plots of required time in one graph   

'''
# Constants
g = 9.81  # gravitation constant
J = 50  # number of iteration
delta_x = 0.02  # x step length
L = 1  # total x length
n_b = 0
H = 0.01
x = np.arange(0, L, delta_x)

delta_t = 0.01      # Time-step
epsilon = delta_t/1000

# time
t1 = 0
t2 = 1
t3 = 4
tend = t3 + epsilon

# initial values for position and speed
u = np.zeros(J)  # speed
average = np.mean(0.002*np.e**((-(x-0.5)**2)/0.05**2))
n = np.ones(J)*H+0.002*np.e**((-(x-0.5)**2)/0.05**2)-np.ones(J)*average  # position
u_p = np.empty(J, float)  # updated speed
u_p[0], u_p[-1] = 0, 0
n_p = np.empty(J, float)  # updated position
n_b = np.zeros(J)

t = 0.0
while t<tend:
    # the number of iteration is total time / time step

    # make the plot at t=0
    if abs(t-t1)<epsilon:
        plt.plot(x, n, label='t=0s')

    # following is the FTCS method
    for i in range(1,J-1):
        # update the position and speed
        u_p[i] = u[i]-0.5*(delta_t/delta_x)*(0.5*(u[i+1])**2+g*n[i+1]-0.5*(u[i-1])**2-g*n[i-1])
        n_p[i] = n[i]-0.5*(delta_t/delta_x)*(u[i+1]*(n[i+1]-n_b[i+1])-u[i-1]*(n[i-1]-n_b[i-1]))

    # update the position and speed at boundary points using forward/backward
    # difference method
    u_p[-1] = u[-1]-(delta_t/delta_x)*(0.5*(u[-1])**2+g*n[-1]-0.5*(u[-2])**2-g*n[-2])
    n_p[-1] = n[-1] - (delta_t / delta_x) * (
                u[-1] * (n[-1] - n_b[-1]) - u[-2] * (n[-2] - n_b[-2]))
    u_p[0] = u[0] - (delta_t / delta_x) * (
                0.5 * (u[1]) ** 2 + g * n[1] - 0.5 * (
        u[0]) ** 2 - g * n[0])
    n_p[0] = n[0] - (delta_t / delta_x) * (
                u[1] * (n[1] - n_b[1]) - u[0] * (
                    n[0] - n_b[0]))

    u, u_p = np.copy(u_p), np.copy(u)
    n, n_p = np.copy(n_p), np.copy(n)
    t += delta_t

    # Make plots at the given times
    if abs(t-t2)<epsilon:
        plt.plot(x, n, label='t=1s')
    if abs(t-t3)<epsilon:
        plt.plot(x, n, label='t=4s')

    # uncomment the following code for animation

    # plt.clf()  # clear the plot
    # plt.plot(x, n, label='wave')  # plot the current sin curve
    # plt.ylim([0, 0.015])
    # plt.plot(x, n_b, label='ground')
    # plt.legend()
    # plt.draw()
    # plt.pause(0.01)  # pause to allow a smooth animation


plt.xlabel("position (m)")
plt.ylabel("amplitude (m)")
plt.title('Wave Behaviour Simulation at Different Time \n'
          'Using FTCS Method')
plt.legend()
plt.show()
