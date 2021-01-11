import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import matplotlib.axes as ax
from time import time
from scipy import signal

# Q1 b
'''
pseudo code:
-Write out all the constants needed for this code, and convert them to correct
unit.
-Write out equations of x/y velocity and x/y position as functions.
-Run a for-loop multiple times, and for each for-loop, run the x/y position 
function, and used the newly obtained x/y position to run the x/y velocity 
function. Then recorder the output for these functions
-Finally, plots using recorded x/y position.

'''
Ms = 1  # unit in solar mass
x = 0.47  # unit in AU
y = 0  # unit in AU
V_x = 0  # unit in AU per year
V_y = 8.17  # unit in AU per year
delta_t = 0.0001  # unit in year
G = 39.5  # unit in AU^3*Ms^-1*yr^-2
a = 0.01  # unit in AU^2

# Q1 c

# equations for updating velocity and position are written as following functions
def velocity_x(x, y, V_x):
    # equation for updating velocity in x direction
    return -G*Ms*x*delta_t/(x**2+y**2)**(3/2)+V_x


def velocity_y(x, y, V_y):
    # equation for updating velocity in y direction
    return -G*Ms*y*delta_t/(x**2+y**2)**(3/2)+V_y


def position_x(x, V_x):
    # equation for updating x position
    return V_x*delta_t + x


def position_y(y, V_y):
    # equation for updating y position
    return V_y*delta_t + y


# these are variables used to stored all position and velocity
recorded_Vx, recorded_Vy, recorded_x, recorded_y = [V_x], [V_y], [x], [y]

for i in range(int(1/delta_t)):
    # number of iteration is 1 year divided by delta t
    # update the velocity first, and then update the position.
    # Store the number after each update

    V_x = velocity_x(x, y, V_x)
    recorded_Vx.append(V_x)

    V_y = velocity_y(x, y, V_y)
    recorded_Vy.append(V_y)

    x = position_x(x, V_x)
    recorded_x.append(x)

    y = position_y(y, V_y)
    recorded_y.append(y)


position_graph = plt.plot(recorded_x, recorded_y, label="Mercury's path")
# plot a graph of position using recorded number for x and y

plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.title('Change in Position of Mercury in a Year Due to \n'
          'Gravitational Force of The Sun')

plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
# To make the axes even by setting their length

plt.legend()
plt.show()


# following are calculation to prove angular momentum conserved
L_diff = Ms*(recorded_x[0]*recorded_Vy[0]-recorded_y[0]*recorded_Vx[0]-
                   (recorded_x[-1]*recorded_Vy[-1]-recorded_y[-1]*recorded_Vx[-1]))

print(L_diff)  # L_diff equal to 0
# The angular momentum difference is 0, this implies the angular momentum
# conserved after the orbiting.


# -----------------------------------------------------------------
# Q1 d

x = 0.47  # unit in AU
y = 0  # unit in AU
V_x = 0  # unit in AU per year
V_y = 8.17  # unit in AU per year


# modify the equations for velocity based on general relativity
def velocity_x_GR(x, y, V_x):
    return -G*Ms*x*delta_t*(1+a/(x**2+y**2))/(x**2+y**2)**(3/2)+V_x


def velocity_y_GR(x, y, V_y):
    return -G*Ms*y*delta_t*(1+a/(x**2+y**2))/(x**2+y**2)**(3/2)+V_y


# these are variables used to stored all position and velocity
recorded_Vx_GR, recorded_Vy_GR, recorded_x_GR, recorded_y_GR = [V_x], [V_y], [x], [y]
recorded_r = [(x**2+y**2)**0.5]  # Store distance bewteen sun and Mercury
for i in range(int(1/delta_t)):

    V_x = velocity_x_GR(x, y, V_x)
    recorded_Vx_GR.append(V_x)

    V_y = velocity_y_GR(x, y, V_y)
    recorded_Vy_GR.append(V_y)

    x = position_x(x, V_x)
    recorded_x_GR.append(x)

    y = position_y(y, V_y)
    recorded_y_GR.append(y)

    recorded_r.append((x**2+y**2)**0.5)

recorded_x_GR, recorded_y_GR = np.array(recorded_x_GR), np.array(recorded_y_GR)
perihelion_position = signal.argrelextrema(np.array(recorded_r), np.less)
# perihelion_position is the index of all local minimum in recorded_r

position_graph_GR = plt.plot(recorded_x_GR, recorded_y_GR, label="Mercury's path")
recorded_x, recorded_y = np.array(recorded_x), np.array(recorded_y)
plt.plot(recorded_x_GR[perihelion_position], recorded_y_GR[perihelion_position],
         'o', label='perihelion')
# plot the perihelion alone with Mercury's position graph, so that it's easier
# to where are the perihelion.

plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.title('Change in Position of Mercury in a Year Due to Gravitational \n'
          'Force of The Sun Predicted by General Relativity')
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
# To make the axes even by setting their length

plt.legend()
plt.show()


# ------------------------------------------------------------------
# Q3

'''
pseudo code:
-Firstly, set the minimum and maximum size of matrix
-Then, run a for-loop, in each iteration, increase the matrix size by one, and 
record the time required for each method
-Finally, write the recorded data in a txt file, and plot the size vs time 
graphs. 
'''
#
starting_N = 2
ending_N = 150
t_method1 = []
# This stored the time required for matrix row multiplication

t_method2 = []
# This stored the time required for matrix multiplication done by numpy.dot


# The following for-loop is from the textbook, it's used to calculate the
# run-time of matrix multiplication if it's being done multiplying row by row.
for N in range(starting_N, ending_N):
    # each iteration has
    C = np.zeros([N, N], float)
    A = np.ones([N, N]) * 3
    B = np.ones([N, N]) * 4

    start = time()
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i, j] += A[i, k] * B[k, j]

    end = time()
    t_method1.append(end - start)
    # record the time required using row multiplication

    start = time()
    C = np.dot(A, B)
    end = time()
    t_method2.append(end - start)
    # record the time required using nump.dot

# Write the run time for each method into a txt file.
f = open('run_time.txt', 'w')
f.write('run time using row multiplication\n')
f.write('matrix size \t run time\n')
for i in range(len(t_method1)):
    f.write('{} \t {}\n'.format(starting_N + i, t_method1[i]))

f.write('\n\n')
f.write('run time using numpy.dot\n')
f.write('matrix size \t run time\n')
for i in range(len(t_method1)):
    f.write('{} \t {}\n'.format(starting_N + i, t_method2[i]))

f.close()

# plot the size vs run-time together for a better comparison.
plot_1 = plt.plot(np.arange(starting_N, ending_N), t_method1,
                  label='row multiplication')
plot_2 = plt.plot(np.arange(starting_N, ending_N), t_method2, label='numpy.dot')
plt.xlabel('Size of Matrix')
plt.ylabel('Time (s)')
plt.title('Increase in run time with respect to increase in matrix size when \n'
          'doing matrix multiplication using for-loop vs using numpy.dot')

plt.show()

'''
From this plot, we can see that the run-time for using row multiplication 
increases like function x^3. It means as matrix size increases, the program will
take cubic time. On the other hand, run-time for using numpy.dot is 
almost 0, this method takes constant time.  
'''
