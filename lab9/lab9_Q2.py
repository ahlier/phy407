import numpy as np
import matplotlib.pyplot as plt

'''
This program is design to investigate the change in electric field and H in a 
square cavity with driving current of different frequency

pseudo-code:
-first define functions for discrete sin/cos transform and inverse discrete 
sin/cos transform
-Write current density function
-Write Crank-Nicolson scheme functions tp update the time on E and H
-write a for-loop, each iteration uses Crank-Nicolson scheme to update the time
and store the value
-make plot using given points of interest.
'''

# defining constant
tao = 0.01  # time interval
Lx = 1  # length of cavity
Ly = 1
J0 = 1
m = 1
n = 1
c = 1
P = 32
T = 20
N = int(T//tao)
D_x = np.pi*c*tao/(2*Lx)
D_y = np.pi*c*tao/(2*Ly)
time = np.linspace(0, T, N)


# The following are fourier transform functions and inverse fourier transform
# functions
def dct(y):
    N = len(y)
    y2 = np.empty(2*N,float)
    y2[:N] = y[:]
    y2[N:] = y[::-1]

    c = np.fft.rfft(y2)
    phi = np.exp(-1j*np.pi*np.arange(N)/(2*N))
    return np.real(phi*c[:N])

def dst(y):
    N = len(y)
    y2 = np.empty(2*N,float)
    y2[0] = y2[N] = 0.0
    y2[1:N] = y[1:]
    y2[:N:-1] = -y[1:]
    a = -np.imag(np.fft.rfft(y2))[:N]
    a[0] = 0.0

    return a

def idst(a):
    N = len(a)
    c = np.empty(N+1,complex)
    c[0] = c[N] = 0.0
    c[1:N] = -1j*a[1:]
    y = np.fft.irfft(c)[:N]
    y[0] = 0.0

    return y

def idst2(b):
    M = b.shape[0]
    N = b.shape[1]
    a = np.empty([M,N],float)
    y = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = idst(b[i,:])
    for j in range(N):
        y[:,j] = idst(a[:,j])

    return y

def idct(a):
    N = len(a)
    c = np.empty(N+1,complex)

    phi = np.exp(1j*np.pi*np.arange(N)/(2*N))
    c[:N] = phi*a
    c[N] = 0.0
    return np.fft.irfft(c)[:N]

def idcst(b):
    # used on function has shape cos()sin()
    M = b.shape[0]
    N = b.shape[1]
    a = np.empty([M,N],float)
    y = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = idct(b[i,:])
    for j in range(N):
        y[:,j] = idst(a[:,j])

    return y

def idsct(b):
    # used on function has shape sin()cos()
    M = b.shape[0]
    N = b.shape[1]
    a = np.empty([M,N],float)
    y = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = idst(b[i,:])
    for j in range(N):
        y[:,j] = idct(a[:,j])

    return y

def dst2(y):
    M = y.shape[0]
    N = y.shape[1]
    a = np.empty([M,N],float)
    b = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = dst(y[i,:])
    for j in range(N):
        b[:,j] = dst(a[:,j])

    return b

def dsct(y):
    # used on function has shape sin()cos()
    M = y.shape[0]
    N = y.shape[1]
    a = np.empty([M,N],float)
    b = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = dst(y[i,:])
    for j in range(N):
        b[:,j] = dct(a[:,j])

    return b

def dcst(y):
    # used on function has shape cos()sin()
    M = y.shape[0]
    N = y.shape[1]
    a = np.empty([M,N],float)
    b = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = dct(y[i,:])
    for j in range(N):
        b[:,j] = dst(a[:,j])

    return b


# The following are Crank-Nicolson scheme for updating time
def update_E_hat(p, q, E_hat, x_hat, y_hat, j_hat):
    return ((1-p**2*D_x**2-q**2*D_y**2)*E_hat+2*q*D_y*x_hat+2*p*D_x*y_hat+tao*j_hat)/(1+p**2*D_x**2+q**2*D_y**2)

def update_x_hat(x_hat, q, updated_E_hat, E_hat):
    return x_hat-q*D_y*(updated_E_hat+E_hat)

def update_y_hat(y_hat, p, updated_E_hat, E_hat):
    return y_hat-p*D_x*(updated_E_hat+E_hat)

# Current density function
def J_z(x, y, t):
    return J0*np.sin(m*np.pi*x/Lx)*np.sin(n*np.pi*y/Ly)*np.sin(w*t)


# This function can be used to create current density matrix at time t
def get_J(t):
    J = np.empty((P, P))
    x = np.linspace(0, Lx, P)
    y = np.linspace(0, Ly, P)
    for i in range(len(y)):
        for k in range(len(x)):
            J[i][k] = J_z(x[k], y[i], t)
    return J


def get_EB():
    '''
    This function will run through the simulation, and return array for Hx, Hy,
    and E at point Hx(0.5, 0), Hy(0, 0.5), and E(0.5, 0.5)
    '''

    # create matrix for E, Hx and Hy
    E = np.zeros((P, P))
    H_x = np.zeros((P, P))
    H_y = np.zeros((P, P))

    # used to record the value of points of interest
    H_x_section = []
    H_y_section = []
    E_section = []
    for i in range(N):
        # each iteration will use discrete sin/cos transform on E, Hx, and Hy
        # then use Crank-Nicolson to update the time, and finally use inverse
        # discrete sin/cos transform to get the E, Hx, and Hy matrix.
        J = get_J(i * tao)
        E_hat = dst2(E)
        x_hat = dsct(H_x)
        y_hat = dcst(H_y)
        J_hat = dst2(J)
        E_hat_f = np.copy(E_hat)
        x_hat_f = np.copy(x_hat)
        y_hat_f = np.copy(y_hat)

        for j in range(P):
            for k in range(P):
                # if statement is to check if the points are at the boundary
                if j != 0 and k != 0 and j != P - 1 and k != P - 1:
                    E_hat_f[j][k] = update_E_hat(k, j, E_hat[j][k], x_hat[j][k],
                                                 y_hat[j][k], J_hat[j][k])

                if j != 0 and j != P - 1:
                    x_hat_f[j][k] = update_x_hat(x_hat[j][k], j, E_hat_f[j][k],
                                                 E_hat[j][k])

                if k != 0 and k != P - 1:
                    y_hat_f[j][k] = update_y_hat(y_hat[j][k], k, E_hat_f[j][k],
                                                 E_hat[j][k])

        E = idst2(E_hat_f)
        H_x = idsct(x_hat_f)
        H_y = idcst(y_hat_f)

        H_x_section.append(H_x[0, int(P // 2)])
        H_y_section.append(H_y[int(P // 2), 0])
        E_section.append(E[int(P // 2), int(P // 2)])
    return H_x_section, H_y_section, E_section


# Q2_c --------------------------------------------------
# make plot
w = 3.75
H_x_section, H_y_section, E_section = get_EB()
plt.plot(time, H_x_section, 'o',label='H_x')
plt.plot(time, H_y_section, '+',label='H_y')
plt.plot(time, E_section, label='E')
plt.xlabel('time (s)')
plt.ylabel('amplitude')
plt.title('Change of E-Field, Hx, and Hy \n'
          'with Respect to Change in Time at {}=3.75'.format('\u03C9'))
plt.legend()
plt.show()



# Q2_d ----------------------------------------------------
# create various value of w, and make plot based on different values of w
w_range = np.linspace(0, 9, 50)
H_x_max, H_y_max, E_max = [], [], []
for i in w_range:
    # get the maximum value of E
    w = i
    H_x_section, H_y_section, E_section = get_EB()
    E_max.append(max(E_section))

plt.plot(w_range, E_max, label='Max of E')
plt.xlabel('Frequency')
plt.ylabel('Electric Field')
plt.title('Change in Maximum Value of Electric Field with \n'
          'Respect to Change in Frequency')
plt.legend()
plt.show()



# Q2_e ----------------------------------------------------
# same as Q2_c
w = np.pi*np.sqrt(2)
H_x_section, H_y_section, E_section = get_EB()
plt.plot(time, H_x_section, 'o',label='H_x')
plt.plot(time, H_y_section, '+',label='H_y')
plt.plot(time, E_section, label='E')
plt.xlabel('time (s)')
plt.ylabel('amplitude')
plt.title('Change of E-Field, Hx, and Hy \n'
          'with Respect to Change in Time at {}=sqrt(2)pi'.format('\u03C9'))
plt.legend()
plt.show()
