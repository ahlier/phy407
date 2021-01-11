import numpy as np
import matplotlib.pyplot as plt

'''
This program will simulate motion of 2 particles under influence of lennard 
jones potential

pseudo code:
- write out equation for acceleration from lennard jones potential
- run a for-loop with N iteration, and for each iteration, apply the Verlet 
algorithm to update position and velocity, and store all the position value
- use stored position values to make plot
'''
m = 1  # mass of particle
sigma = 1  # given
epsilon = 1  # given
dt = 0.01  # time interval
N = 100  # total step taken
point_set = [[[4, 4], [5.2, 4]], [[4.5, 4], [5.2, 4]], [[2, 3], [3.5, 4.4]]]
# each element in point_set represents one scenario, each element contains
# two points
def acceleration(x, y):
    '''
    This function will calculate acceleration of a particle based on lennard
    jones potential equation

    :param x: x position of particle
    :param y: y position of particle
    :return: x y component of acceleration
    '''
    a = 4*epsilon*(12*sigma**12/(x**2+y**2)**(13/2)-6*sigma**6/(x**2+y**2)**(7/2))/m
    a_x = a*x/np.sqrt(x**2+y**2)
    a_y = a*y/np.sqrt(x**2+y**2)
    return a_x, a_y


def verlet(dt, x, y, v_x, v_y):
    '''
    This is the verlet algorithm, it is used to update position and velocity
    after a small time period dt

    :param dt: time interval
    :param x: x position
    :param y: y position
    :param v_x: x component of velocity
    :param v_y: y component of velocity
    :return: updated x,y position and x,y component of velocity
    '''
    new_x = x+dt*v_x
    new_y = y+dt*v_y

    a_x, a_y = acceleration(new_x, new_y)
    k_x, k_y = dt*a_x, dt*a_y
    v_new_x = v_x+k_x
    v_new_y = v_y+k_y
    return new_x, new_y, v_new_x, v_new_y


for p in point_set:
    # each iteration will run a scenario of 2 particles, p contains position of
    # these two particle
    x = p[0][0]-p[1][0]  # x separation
    y = p[0][1]-p[1][1]  # y separation

    p1_x, p1_y = [p[0][0]], [p[0][1]]  # x,y position of first particle
    p2_x, p2_y = [p[1][0]], [p[1][1]]

    a_x, a_y = acceleration(x, y)  # x,y component of acceleration
    v_x = 0.5*dt*a_x  # x velocity in first dt/2 sec
    v_y = 0.5*dt*a_y  # y velocity in first dt/2 sec

    for i in range(N):
        # each iteration will use verlet algorithm to update the position and
        # velocity, and also store the updated position

        x, y, v_x, v_y = verlet(dt, x, y, v_x, v_y)

        delta_x = x-(p1_x[-1]-p2_x[-1])  # change in x direction
        delta_y = y-(p1_y[-1]-p2_y[-1])  # change in y direction

        # The following update the particle position
        p1_x_new = 0.5 * delta_x + p1_x[-1]
        p2_x_new = p2_x[-1] - 0.5 * delta_x

        p1_y_new = 0.5 * delta_y + p1_y[-1]
        p2_y_new = p2_y[-1] - 0.5 * delta_y

        p1_x.append(p1_x_new)
        p1_y.append(p1_y_new)
        p2_x.append(p2_x_new)
        p2_y.append(p2_y_new)

    plt.plot(p1_x, p1_y, '.', label='point 1')
    plt.plot(p2_x, p2_y, '.', label='point 2')
    plt.title('Motion of Two particles Described by lennard Jones Potential')
    plt.legend()

    plt.show()



