import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
import time

"""
Question 1 - Part D
Purpose - Code the Crank-Nicolson for the wave equation in the time-dependent Schrodinger equation in Double Well

# Pseudo code
# Set constants
#   Electron mass
#   Total Length
#   Initial localization
#   Sigma
#   K
#   Planck's constant
#   THe number of time steps
#   Time step
#   The number of positions
#   The length between each position
# Define the wave function - Eq 20 from instruction
# Create an array to hold values
#   positions
#   wave functions
# Normalize the initial wave function
# Define potential and constants in Hamiltonian (A, B)
#   V = 0
# Create the discretized Hamiltonian
# Create the identity matrix
# Create the matrices L and R as defined in equation 17
# Compute the wave functions at each time step
# Compute average position at each time step
# Plot the graphs
#   Average position
"""

# Constants
m = 9.109e-31 # Electron mass, Unit: kg
L = 1e-8
x0 = L/3
x1 = L/4
sigma = L/25
k = 500/L

h_bar = 6.626 * 10 ** (-34)     # Planck constant
N = 6000                        # Number of time steps
P = 1024                        # Number of spatial segments
a = L / P                       # Distance between points
t = 10 ** (-18)                 # time step in seconds
V0 = 6 * 10**(-17)              # Initial Potential
ftsz = 24

def psi_0(x):
    return np.exp(-(x - x0) ** 2 / (4 * sigma ** 2) + 1j * k * x)


# Define space distribution and psi
#x_points = np.linspace(-L/2+a, (P-1)*a-L/2, P)
x_points = np.linspace(a-L/2, (P-1)*a-L/2, P-1)
psi = np.array(list(map(psi_0, x_points)), complex)
normalization = np.sqrt(np.matmul(np.conj(psi), psi))
psi_normalized = psi / normalization

# Define potential
V_d = V0*((x_points)**2/(x1)**2 - 1)**2 # Potential for Double well

# Define the discretized Hamiltonian
A = -1*(h_bar**2)/(2*m*(a**2))
B = V_d + (-2*A)
H_D = np.zeros([P-1, P-1], complex)

for i in range(len(H_D)):
    if i == 0:
        H_D[i][i] = B[i]
        H_D[i][i+1] = A

    elif i == (len(H_D)-1):
        H_D[i][i-1] = A
        H_D[i][i] = B[i]

    else:
        H_D[i][i-1] = A
        H_D[i][i] = B[i]
        H_D[i][i+1] = A

# Define P-1 identity matrix
I = np.zeros([P-1, P-1], complex)
for i in range(len(I)):
    I[i][i] = 1+0j

# Compute the matrix L and R
L = np.zeros([P-1, P-1], complex)
R = np.zeros([P-1, P-1], complex)
for i in range(len(L)):
    for j in range(len(L)):
        L[i][j] = I[i][j] + 1j * t * H_D[i][j] / (2*h_bar)
        R[i][j] = I[i][j] - 1j * t * H_D[i][j] / (2*h_bar)

# Time evolution
solution = []
solution.append(psi_normalized)

# Calculate the wave equation at each time step
for i in range(N):
    v = np.matmul(R, solution[-1])
    psi_new = solve(L, v)
    normalization = np.sqrt(np.matmul(np.conj(psi_new), psi_new))
    psi_new_normalized = psi_new / normalization
    solution.append(psi_new_normalized)

t_array = []
position_array = []
for i in range(N):
    sum_position = np.matmul(x_points, abs(solution[i]) ** 2)           #| psi|^2 * x
    average_position = sum_position / len(x_points)                     # Average Position

    t_array.append(t * i)
    position_array.append(average_position)

# Plot of average position over time
plt.figure(1)
plt.plot(t_array, position_array, 'b-',label="Position")
plt.xlabel("Time (in second)", fontsize=ftsz)
plt.ylabel("Position (in meter)", fontsize=ftsz)
plt.title("Position of the wave function with time, <X>(t), Double well", fontsize=ftsz)
plt.legend(fontsize=ftsz)

"""
Question 1 - Part d
Purpose - Create an animation of the wave function inside Double well

# Pseudo code
# (Follow the structure of the sample code for animation by Professor)
# Create Figure
# Create axis
# Define functions
#   Init
#   animate
# Create animation by using FuncAnimation (interval 1 ensures the fastest animation)
"""

from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-deep')

fig = plt.figure(2)
ax = plt.axes(xlim=(x_points[0], x_points[-1]), ylim=(0, 0.012), xlabel="position (in meter)", ylabel="$|\psi|^2$", title="Probability Density $|\psi|^2$ of the wave function in a Double well")
line, = ax.plot([], [], lw=3)

def init():
    line.set_data([],[])
    return line,

def animate(i):
    x = x_points
    y = abs(solution[i])**2
    line.set_data(x, y)
    return line,

anim = FuncAnimation(fig, animate, init_func=init, frames=N, interval=1, blit=True)
#anim.save('Probability_Density_DoubleWell.gif')
plt.show()