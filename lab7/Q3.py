import matplotlib.pyplot as plt
import scipy.constants as pc
import numpy as np

'''
This program can calculate the eigen energy of hydrogen atom with given n and l,
it can also plot the R and R^2 function

pseudo code:
-Write out the potential function
-write out a function f that is able to return the derivative of the pair of ODEs
-Use 4th order Runge-Kutta method to get R value in a given x range
-Use secant method to get the eigen energy and R function
-normalized R function using Simpson's rule
-make plots
'''

e = pc.e  # electron charge
a = 5e-11  # in m
hbar = pc.hbar
m = pc.electron_mass  # electron mass in kg
h = 0.002*a  # step size
x_0 = h  # initial value of x
x_f = 20*a  # final value of x
l = 0  # order
N = int((x_f-x_0)//h)  # number of steps


def psi(E):
    """
    This function uses 4th order Runge-Kutta to calculate value of R for a given
    range of x.
    :param E: energy
    :return: the R function
    """
    def f(r, x):
        '''
        :param r: it contains value of R and S
        :param x: position
        :return: derivative of R and S
        '''

        def V(x):
            # potential function
            return -(e**2)/(4*pc.pi*pc.epsilon_0*x)

        R = r[0]
        S = r[1]
        f_R = S
        f_S = (l*(l+1)*R+2*m*x**2/hbar**2*(V(x)-E)*R-2*x*S)/x**2
        return np.array([f_R, f_S], float)

    r = np.array([0, 1] ,float) # r is the value of R and S
    R_function = []
    for x in np.arange(x_0, x_f, h):
        # 4th order Runge-Kutta method
        R_function.append(r[0])
        k1 = h * f(r, x)
        k2 = h * f(r + 0.5 * k1, x + 0.5 * h)
        k3 = h * f(r + 0.5 * k2, x + 0.5 * h)
        k4 = h * f(r + k3, x + h)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return np.array(R_function, float)


def secant(E1, E2):
    '''
    :param E1: energy 1
    :param E2: energy 2
    :return: eigen energy and normalized R function
    '''
    target_accuracy = e / 100000  # Used to check if the function converges
    R = psi(E1)
    psi2 = R[- 1]
    while abs(E1 - E2) > target_accuracy:
        # secant method
        R = psi(E2)
        psi1, psi2 = psi2, R[- 1]
        E1, E2 = E2, E2 - psi2 * (E2 - E1) / (psi2 - psi1)
    # Normalizing using Simpson's rule
    R_square = R * R
    integral = h / 3 *(R_square[0] + R_square[N//2 - 1] + \
            4 * sum(R_square[1 : N//2 : 2]) + 2 * sum(R_square[0 : N//2 + 1 : 2]) )

    return E2 / e, R/np.sqrt(2*integral)


a_0 = 0.0529*10**(-9)  # in m
def r_1(x):
    # theoretical function for R when n=1
    return np.e**(-x/a_0)/a_0

def r_2(x):
    # theoretical function for R when n=2
    return (2-x/a_0)*np.e**(-x/(2*a_0))/(10**10*np.sqrt(8)*a_0**1.5)

E1, psi1 = secant(-15*e, -13*e)
E2, psi2 = secant(-15*e/4, -13*e/4)
print('E_1 = ', E1, 'eV')
print('E_2 = ', E2, 'eV')

# graphing
xpoints = np.arange(x_0, x_f, h)

# plt.plot(xpoints, psi1, label='R of n=1, l=0')
# plt.plot(xpoints, psi2, label='R of n=2, l=0')
l = 1

E, psi3 = secant(-15*e/4, -13*e/4)
plt.plot(xpoints, psi3, label='R of n=2, l=1')
print(print('E = ', E, 'eV'))

plt.xlabel('x (m)')
plt.title('Normalized R Function of Hydrogen Atom')
plt.ylabel('R')
plt.legend()
plt.show()



plt.plot(xpoints, psi1**2, label='R^2 of n=1, l=0')
plt.plot(xpoints, psi2**2, label='R^2 of n=2, l=0')
plt.plot(xpoints, psi3**2, label='R^2 of n=2, l=1')
plt.xlabel('x (m)')
plt.title('Normalized R^2 Function of Hydrogen Atom')
plt.ylabel('R^2')
plt.legend()
plt.show()


plt.plot(xpoints, psi1**2, label='Numerical R')
plt.plot(xpoints, [r_1(i) for i in xpoints], label='Theoretical R')
plt.xlabel('x (m)')
plt.title('Comparison Between Numerical and \n'
          'Theoretical Value of R when n=1')
plt.ylabel('R')
plt.legend()
plt.show()


plt.plot(xpoints, psi2, label='Numerical R')
plt.plot(xpoints, [r_2(i) for i in xpoints], label='Theoretical R')
plt.xlabel('x (m)')
plt.title('Comparison Between Numerical and \n'
          'Theoretical Value of R when n=2')
plt.ylabel('R')
plt.legend()
plt.show()
