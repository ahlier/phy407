import numpy as np
import matplotlib.pyplot as plt


'''
This program is designed to investigated the behaviour of wave function in an
infinite potential well with non-constant potential across the well

pseudo code:
- create function to calculate H value with given m and n, and psi function
- using 2 for-loops to find the H function with given size
- Find the eigenvalue of H function and use it to calculate probability density
- Calculate the integral using trapezoidal rule
'''
# Q2_b
h_bar = 6.582119514 * 10 ** -16 # in J * s
L = 5 * 10 ** -10 # in m
M = 9.1094 * 10 ** -31 # in kg
a = 10 # in eV

def find_Hmn(m, n):
    # The H function is given, the output is base on m and n, which are the
    # position of element
    if (m+n)%2 == 1:
        return - 8 * a / (np.pi ** 2)  * m * n / (m ** 2 - n ** 2) ** 2

    elif m == n:
        return (np.pi ** 2 * h_bar ** 2 * m ** 2) / (2 * M * L ** 2)*1.6*10**(-19)+ a / 2
    else:
        return 0

# Q2_c
def get_H(size):
    # general H function with given size
    H = np.zeros((size, size))
    for m in range(size):
        for n in range(size):
            H[m][n] = find_Hmn(m+1, n+1)
    return H

# Find the lowest eigenvalue of H of size 10
size = 10
H_matrix = get_H(size)
ev = np.linalg.eigvalsh(H_matrix)  # get the eigenvalue of H
ev.sort()
print('The ground state energy is the smallest eigenvalue {}'.format(ev[0]))

# Q2_d
# Find the ten lowest eigenvalue of H of size 100
size = 100
H_matrix = get_H(size)
ev = np.linalg.eigvalsh(H_matrix)  # get the eigenvalue of H
ev.sort()
print('The first tem energy eigenvalues are {}'.format(ev[:10]))


# Q2_e

e_value, e_vector = np.linalg.eigh(H_matrix)


def psi(n, x):
    # wave function is given
    amplitude = 0
    for m in range(size):
        amplitude += (2 / L)**(1/2) * e_vector[n][m] * np.sin(np.pi * x * (m+1) / L)
    return amplitude


# plot three lowest energy states
x = np.linspace(0, L, 100)
v = [i*a/L for i in x]

# probability density
psi_0_squared = [(psi(0, i))**2 for i in x]
psi_1_squared = [(psi(1, i))**2 for i in x]
psi_2_squared = [(psi(2, i))**2 for i in x]

plt.plot(x, psi_0_squared, label='psi0')
plt.plot(x, psi_1_squared, label='psi1')
plt.plot(x, psi_2_squared, label='psi2')
plt.xlabel('position (m)')
plt.ylabel('probability density')
plt.title('Probability Density of First Three States in an \n'
          'Infinite Square Well')
plt.legend()
plt.show()


# The following are probability density function first three states
def psi_0_square(x):
    return psi(0, x)**2


def psi_1_square(x):
    return psi(1, x)**2


def psi_2_square(x):
    return psi(2, x)**2


def trap(psi_n_square, n):
    # The integral is done using trapezoidal rule, range is 0 to L
    h = 1 / float(n)
    intgr = 0.5 * h * (psi_n_square(0) + psi_n_square(L))
    for i in range(1, int(n)):
        intgr += h * psi_n_square(i * h)
    return intgr

# Calculating the integral from 0 to L
print(trap(psi_0_square, 100)/np.sqrt(L))
print(trap(psi_1_square, 100)/np.sqrt(L))
print(trap(psi_2_square, 100)/np.sqrt(L))
