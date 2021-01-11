import numpy as np
import matplotlib.pyplot as plt


'''
This program is similar to Q2, it is also used to simulate behaviour of a wave,
but it uses Two-Step Lax-Wendroff scheme instead of FTCS.

pseudo code:
-First set the values constants and variables
-Then set the function for force
-Run a loop, each iteration will use Two-Step Lax-Wendroff method to update the 
position and speed of the wave, also update the boundary points using 
forward/backward difference method.
-plot the graph in each iteration to create an animation
-make the plots of required time in one graph  
'''
# Constants, same as Q2
g = 9.81
J = 50
delta_x = 0.02
L = 1
H = 0.01
x = np.arange(0, L, delta_x)

delta_t = 0.01
epsilon = delta_t / 1000

# following are time of interest
t1 = 0
t2 = 1
t3 = 4
tend = t3 + epsilon

u = np.zeros(J)
average = np.mean(0.002 * np.e ** ((-(x - 0.5) ** 2) / 0.05 ** 2))
n = np.ones(J) * H + 0.002 * np.e ** ((-(x - 0.5) ** 2) / 0.05 ** 2) - np.ones(
    J) * average
u_p = np.empty(J, float)
u_p[0], u_p[-1] = 0, 0
n_p = np.empty(J, float)
n_b = np.zeros(J)

t = 0.0
while t < tend:
    if abs(t-t1)<epsilon:
        plt.plot(x ,n, label='t=0s')

    def F_u(u, n):
        # this returns x-direction of the force
        return 0.5*u**2+g*n

    def F_n(u, n, n_b):
        # this returns y-direction of the force
        return (n-n_b)*u

    for i in range(1, J - 1):
        # The following uses Two-Step Lax-Wendroff method to calculate the
        # values at half point of position and speed, and use them to update
        # position and speed
        u_posi_half = 0.5 * (u[i + 1] + u[i]) - (delta_t / delta_x) * 0.5 * \
                      (F_u(u[i+1], n[i+1])-F_u(u[i], n[i]))
        n_posi_half = 0.5 * (n[i + 1] + n[i]) - (delta_t / delta_x) * 0.5 * \
                      (F_n(u[i+1], n[i+1], n_b[i+1])-F_n(u[i], n[i], n_b[i]))
        u_neg_half = 0.5 * (u[i] + u[i-1]) - (delta_t / delta_x) * 0.5 * \
                     (F_u(u[i], n[i])-F_u(u[i-1], n[i-1]))
        n_neg_half = 0.5 * (n[i] + n[i-1]) - (delta_t / delta_x) * 0.5 * (
                    F_n(u[i], n[i], n_b[i + 1]) - F_n(u[i-1], n[i-1],
                                                              n_b[i-1]))

        u_p[i] = u[i] - (delta_t / delta_x) * (F_u(u_posi_half, n_posi_half) -
                                               F_u(u_neg_half, n_neg_half))
        n_p[i] = n[i] - (delta_t / delta_x) * (F_n(u_posi_half, n_posi_half, n_b[i])-
                                               F_n(u_neg_half, n_neg_half, n_b[i]))

    # update the position and speed at boundary points using forward/backward
    # difference method
    u_p[-1] = u[-1] - (delta_t / delta_x) * (
                0.5 * (u[-1]) ** 2 + g * n[-1] - 0.5 * (u[-2]) ** 2 - g * n[-2])
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

    # uncomment the following for animation

    # plt.clf()  # clear the plot
    # plt.plot(x, n, label='wave')  # plot the current sin curve
    # plt.ylim([0, 0.015])
    # plt.plot(x, n_b, label='ground')
    # plt.legend()
    # plt.draw()
    # plt.pause(0.01)  # pause to allow a smooth animation

plt.xlabel("position (m)")
plt.ylabel("amplitude (m)")
plt.title('Wave Behaviour Simulation at Different Time Using \n'
          'Two-Step Lax-Wendroff method')
plt.legend()
plt.show()




#-----------------------------------------------------------------------
# Q3 part b
'''
part b uses code from part a, but value of some variables are changed, such as 
the ground position. Also, the boundary points are set to be rigid
'''

# Constants
J = 150
delta_x = 1/J
L = 1
x = np.arange(0, L, delta_x)

delta_t = 0.001  # Time-step
epsilon = delta_t / 1000

t1 = 0
t2 = 1
t3 = 4
t4 = 5  # because at least 5000 iterations are required
tend = t4 + epsilon

# create initial value array for position of wave and ground, and wave speed
u = np.zeros(J)
average = np.mean(0.002 * np.e ** ((-(x - 0.5) ** 2) / 0.05 ** 2))
n = np.ones(J) * H + 0.002 * np.e ** ((-(x - 0.5) ** 2) / 0.05 ** 2) - np.ones(
    J) * average
u_p = np.empty(J, float)
u_p[0], u_p[-1] = 0, 0
n_p = np.empty(J, float)
n_bs = H-0.0004
n_b = n_bs*0.5*(1+np.tanh(180*(x-0.5)/(np.pi*8*np.pi)))  # position of ground

t = 0.0
# The boundary is rigid, so the speed at boundary points are kept the same
boundary_point = [np.copy(n[0]), np.copy(n[-1])]
while t < tend:
    # This loop is the same as the one from part a
    if abs(t-t1)<epsilon:
        plt.plot(x, n, label='t=0s')

    def F_u(u, n):
        return 0.5*u**2+g*n

    def F_n(u, n, n_b):
        return (n-n_b)*u

    n_p[-1] = boundary_point[0]
    n_p[0] = boundary_point[1]

    for i in range(1, J - 1):
        u_posi_half = 0.5 * (u[i + 1] + u[i]) - (delta_t / delta_x) * 0.5 * (F_u(u[i+1], n[i+1])-F_u(u[i], n[i]))
        n_posi_half = 0.5 * (n[i + 1] + n[i]) - (delta_t / delta_x) * 0.5 * (F_n(u[i+1], n[i+1], n_b[i+1])-F_n(u[i], n[i], n_b[i]))
        u_neg_half = 0.5 * (u[i] + u[i-1]) - (delta_t / delta_x) * 0.5 * (F_u(u[i], n[i])-F_u(u[i-1], n[i-1]))
        n_neg_half = 0.5 * (n[i] + n[i-1]) - (delta_t / delta_x) * 0.5 * (
                    F_n(u[i], n[i], n_b[i + 1]) - F_n(u[i-1], n[i-1],
                                                              n_b[i-1]))

        u_p[i] = u[i] - (delta_t / delta_x) * (F_u(u_posi_half, n_posi_half) - F_u(u_neg_half, n_neg_half))
        n_p[i] = n[i] - (delta_t / delta_x) * (F_n(u_posi_half, n_posi_half, n_b[i])-F_n(u_neg_half, n_neg_half, n_b[i]))
    u_p[-1] = 0
    u_p[0] = 0


    u, u_p = np.copy(u_p), np.copy(u)
    n, n_p = np.copy(n_p), np.copy(n)
    t += delta_t

    # Make plots at the given times
    if abs(t-t2)<epsilon:
        plt.plot(x, n, label='t=1s')
    if abs(t-t3)<epsilon:
        plt.plot(x, n, label='t=4s')

    # Uncomment the following for animation

    # plt.clf()  # clear the plot
    # plt.plot(x, n, label='wave')  # plot the current sin curve
    # plt.ylim([0, 0.015])
    # plt.plot(x, n_b, label='ground')
    # plt.legend()
    # plt.draw()
    # plt.pause(0.01)  # pause to allow a smooth animation

plt.plot(x, n_b, label='ground')
plt.xlabel("position (m)")
plt.ylabel("amplitude (m)")
plt.title('Wave Behaviour Simulation at Different Time Using \n'
          'Two-Step Lax-Wendroff method on Uneven Ground')
plt.legend()
plt.show()
