from numpy import copy, empty
import numpy as np
from time import time
import matplotlib.pyplot as plt

'''
This program is designed to investigate the difference between methods; gaussian
elimination, partial pivoting, and LU decomposition when solving matrix from 
perspective of run time and error

pseudo code:
- Create function for methods; gaussian elimination, partial pivoting, and LU 
decomposition.
- generate random matrix of different sizes and test their run time and error
for each method.
- plot the size vs run time and size vs error graphs
'''


# GaussElim function is given
def GaussElim(A_in, v_in):
   """Implement Gaussian Elimination. This should be non-destructive for input
   arrays, so we will copy A and v to
   temporary variables
   IN:
   A_in, the matrix to pivot and triangularize
   v_in, the RHS vector
   OUT:
   x, the vector solution of A_in x = v_in """
   # copy A and v to temporary variables using copy command
   A = copy(A_in)
   v = copy(v_in)
   N = len(v)
   for m in range(N):
       # Divide by the diagonal element
       div = A[m, m]
       if div != 0:
           A[m, :] /= div
           v[m] /= div
        # Now subtract from the lower rows
       for i in range(m + 1, N):
           mult = A[i, m]
           A[i, :] -= mult * A[m, :]
           v[i] -= mult * v[m]
   # Backsubstitution
   # create an array of the same type as the input array
   x = empty(N, dtype=v.dtype)
   for m in range(N-1, -1, -1):
       x[m] = v[m]
       for i in range(m+1, N):
           x[m] -= A[m, i]*x[i]
   return x


def partial_pivoting(A_in, v_in):
    '''
    This program will solve the matrix using partial pivoting method
    :param A_in: input matrix
    :param v_in: output matrix
    :return: solution
    '''
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)

    for m in range(N):

        for i in range(m + 1, N):
            # swap row if the element is larger
            if A[m][m] < A[i][m]:
                v[m], v[i] = copy(v[i]), copy(v[m])
                A[m, :], A[i, :] = copy(A[i, :]), copy(A[m, :])


        # Following is Backsubstitution code from textbook
        div = A[m, m]
        v[m] /= div
        A[m, :] /= div

        for i in range(m + 1, N):
            mult = A[i, m]
            A[i, :] -= mult * A[m, :]
            v[i] -= mult * v[m]

    x = empty(N, float)
    for m in range(N - 1, -1, -1):
        x[m] = v[m]
        for i in range(m + 1, N):
            x[m] -= A[m, i] * x[i]

    return x

def LU_decomposition(A_in, v_in):
    # np.linalg.solve uses LU decomposition
    return np.linalg.solve(A_in, v_in)


size = 200  # max size of randomly generated matrix
gaussian_time, partial_time, LU_time = [], [], []  # time required
gaussian_error, partial_error, LU_error =[], [], []  # error

for i in range(5, size):
    # each iteration will compute the run time and error for size i matrix

    # generate a random matrix of size i, value between 0 and 1
    ran_A = np.random.random((i, i))
    ran_v = np.random.random(i)

    # calculating run time
    start_time_g = time()
    gaussian = GaussElim(ran_A, ran_v)
    end_time_g = time()
    gaussian_time.append(end_time_g-start_time_g)

    start_time_p = time()
    partial = partial_pivoting(ran_A, ran_v)
    end_time_p = time()
    partial_time.append(end_time_p-start_time_p)

    start_time_LU = time()
    LU = LU_decomposition(ran_A, ran_v)
    end_time_LU = time()
    LU_time.append(end_time_LU-start_time_LU)

    # calculating log of error
    gaussian_error.append(np.log(np.mean(abs(np.dot(ran_A, gaussian)-ran_v))))
    partial_error.append(np.log(np.mean(abs(np.dot(ran_A, partial) - ran_v))))
    LU_error.append(np.log(np.mean(abs(np.dot(ran_A, LU) - ran_v))))

# graphing time vs size
plt.plot(np.arange(5, size), gaussian_time, label='gaussian')
plt.plot(np.arange(5, size), partial_time, label='partial')
plt.plot(np.arange(5, size), LU_time, label='LU')
plt.xlabel('size')
plt.ylabel('run time(s)')
plt.title('Increase in Run Time with Respect to Increase in Matrix Size \n'
          'Using Different Methods')
plt.legend()
plt.show()


# graphing log(error) vs log(size)
plt.plot(np.log(np.arange(5, size)), gaussian_error, label='gaussian')
plt.plot(np.log(np.arange(5, size)), partial_error, label='partial')
plt.plot(np.log(np.arange(5, size)), LU_error, label='LU')
plt.xlabel('log(size)')
plt.ylabel('log(error)')
plt.title('Change in Error with Respect to Change in Matrix Size\n'
          'Using Different Methods')
plt.legend()
plt.show()


# Q1_c
R1, R3, R5, R2, R4, R6 = 1000, 1000, 1000, 2000, 2000, 2000 # resistance in ohms
C1, C2 = 10 ** -6, 0.5 * 10 ** -6 # capacitance in farads
xp = 3 # in volts
w = 1000 # in hertz

# The three equations given, but in matrix form, solution will be x value
A = np.array([ [ 1 / R1 + 1 / R4 + 1j * w * C1, - 1j * w * C1, 0 ],
            [ -1j * C1, 1 / R2 + 1 / R5 + 1j * w * C1 + 1j * w * C2, - 1j * w * C2 ],
            [ 0, -1j * w *C2, 1 / R3 + 1/ R6 + 1j * w * C2 ]], complex)

v = np.array([ xp / R1, xp / R2, xp/ R3 ], complex)

x = np.linalg.solve(A, v)
amplitude, phase = [], []
for i in range(len(x)):
    # calculating amplitude and phase from solution
    amplitude.append(abs(x[i]))
    phase.append(np.angle(x[i]))
    print('The amplitude is {}, phase at t=0 is {}'.format(amplitude[-1], phase[-1]*180/np.pi))

# graphing voltage vs time
x = np.linspace(0, 4*np.pi/w, 500)
for i in range(len(amplitude)):

    plt.plot(x, amplitude[i]*np.cos(w*x-phase[i]), label='V{}'.format(i+1))
plt.xlabel('time (s)')
plt.ylabel('Voltage (V)')
plt.title('Change in Voltage with Respect to Change in Time')
plt.legend()
plt.show()
