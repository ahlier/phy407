from numpy import ones, copy, cos, tan, pi, linspace
import numpy as np
import math
import matplotlib.pyplot as plt


# gaussxw and gaussxwab functions are taken from textbook
def gaussxw(N):
    a = linspace(3, 4 * N - 1, N) / (4 * N + 2)
    x = cos(pi * a + 1 / (8 * N * N * tan(a)))

    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = ones(N, float)
        p1 = copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = max(abs(dx))

    w = 2 * (N + 1) * (N + 1) / (N * N * (1 - x * x) * dp * dp)
    return x, w


def gaussxwab(N, a, b):
    x, w = gaussxw(N)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


'''
This program shows the wave function of a particle at different energy levels.
It also show the relation between particle energy and energy level, uncertainty
in position and momentum.


pseudo code:
- To create function H, we need to use recursion, in this recursion, we first 
need to check if n is equal to 0 or 1, if not, then we return a combination of 
H functions with smaller n as input.
- Then we need to create a wave function (the formula is give). 
- To solve Q2a using 2 for-loop, first one being the number of graph we need to 
draw, then second one is the x-range of the graph.
- To solve Q2b, we just need to change n and x-range for wave function
- To solve Q2c, we need to create functions for derivative of psi, mean square 
position/momentum. To find the improper integral in mean square position/momentum,
use change in variable, convert +/- infinity to +/- 1. 
- After we have all the functions, we just need to use a for-loop to run through
n = 0 to n = 15 and make graph of n vs E and uncertainty in position vs momentum.

'''


def H(x, n):
    '''
    recursion is used to find Hermite polynomial

    :param x: position of interest
    :param n: energy level
    :return: Hermite polynomial at given x and n.
    '''
    if n == 1:
        return 2 * x
    elif n == 0:
        return x ** 0
    else:
        return 2 * x * H(x, n - 1) - 2 * (n - 1) * H(x, n - 2)


def psi(x, n):
    '''
    This a the given wave function, a and b are used to avoid overflow.
    :param x: position of interest
    :param n: energy level
    :return: wave amplitude at given x and n
    '''
    a = np.sqrt(2 ** n)
    b = np.e ** (- x ** 2 / 2) * H(x, n) / np.sqrt(math.factorial(n) *
                                                   np.sqrt(np.pi))
    return b / a


# Q2_a-----------------------------------------------------------
# create 4 wave function graphs with n = 0 to 3, and x ranging from -4 to 4
x = np.arange(-4, 4, 0.01)
n = np.arange(0, 4)

for i in n:
    plt.plot(x, psi(x, i), label='n={}'.format(i))
    plt.xlabel('position')
    plt.ylabel('amplitude')
plt.title('Wave Function of a Particle with Different Energy Level')
plt.legend()
plt.show()

# Q2_b------------------------------------------------------------
# This create a wave function with n=30 and x ranging from -10 to 10
N = 1000  # number of points on graph
x_range = np.linspace(-10, 10, N)
H_30 = psi(x_range, 30)
plt.plot(x_range, H_30, label='n=30')
plt.xlabel('position')
plt.ylabel('amplitude')
plt.title('Wave Function of a Particle with Energy Level Equal to 30')
plt.legend()
plt.show()

# Q2_c -----------------------------------------------------
N = 100  # number of sample points using Gaussian quadrature


def mean_square_x(n):
    '''
    This function finds the mean square position with given energy level, the
    improper integral is done by change of variable, make x = z/(1-z^2), thus
    the integral goes from -1 to 1.
    :param n: energy level
    :return: mean square position
    '''

    def z_to_x(z):
        # change of variable
        return z / (1 - z ** 2)

    def integrand(z):
        # (1+z^2)/(1-z^2)^2 is added because of the change of variable
        return (z_to_x(z) ** 2 * (1 + z ** 2) * np.abs(
            psi(z_to_x(z), n)) ** 2) / (1 - z ** 2) ** 2

    integral = 0.0  # value of integral
    x, w = gaussxwab(N, -1, 1)  # x is position, w is weight

    for k in range(N):
        # calculate the integral over each sample point
        integral += w[k] * integrand(x[k])

    return integral


def psi_derivative(x, n):
    '''
    This function calculate the derivative of wave function, n has to be greater
    than 0, because H is not defined when n=-1. a and b are used to avoid
    overflow.
    :param x: position of interest
    :param n: energy level
    :return: derivative of psi at given x, n
    '''
    a = np.sqrt(2 ** n)
    b = np.e ** (- x ** 2 / 2) * (-x * H(x, n) + 2 * n * H(x, n - 1)) / \
        np.sqrt(math.factorial(n) * np.sqrt(np.pi))
    return b / a


def mean_square_p(n):
    '''
    This function is similar to mean_square_x(), this one return the mean square
    momentum at energy level n
    :param n: energy level
    :return: mean square momentum
    '''

    def z_to_x(z):
        # change of variable
        return z / (1 - z ** 2)

    def integrand(z):
        # This integrand is done using change of variable

        return ((1 + z ** 2) * np.abs(psi_derivative(z_to_x(z), n)) ** 2) / \
               (1 - z ** 2) ** 2

    integral = 0.0  # value of integral
    x, w = gaussxwab(N, -1, 1)  # x is position, w is weight

    for k in range(N):
        # calculate the integral over all sample points
        integral += w[k] * integrand(x[k])

    return integral


n_range = np.arange(1, 16)
# energy level of interest is from 1 to 15, n=0 is ignored because psi
# derivative cannot have n=0

mean_sq_x = []  # record mean square position for all n
rt_mean_sq_x = []  # record root mean square position for all n
mean_sq_p = []  # record mean square momentum for all n
rt_mean_sq_p = []  # record root mean square momentum for all n
E = []  # record energy for all n
for i in n_range:
    mean_sq_x.append(mean_square_x(i))
    mean_sq_p.append(mean_square_p(i))
    rt_mean_sq_x.append(np.sqrt(mean_sq_x[-1]))
    rt_mean_sq_p.append(np.sqrt(mean_sq_p[-1]))
    E.append(0.5 * (mean_sq_x[-1] + mean_sq_p[-1]))

# root mean square position/momentum are also uncertainty in position/momentum
# so plot of uncertainty in position vs momentum is:
plt.plot(rt_mean_sq_x, rt_mean_sq_p)
plt.xlabel('position uncertainty')
plt.ylabel('momentum uncertainty')
plt.title('Relation Between Position Uncertainty and Momentum Uncertainty')
plt.show()

# plotting graph of energy energy vs n
plt.plot(np.arange(1, 16), E)
plt.xlabel('energy level')
plt.ylabel('energy (J)')
plt.title(
    'Change in Energy of a Particle with Respect to Change in Energy Level')
plt.show()



