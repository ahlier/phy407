import numpy as np
import matplotlib.pyplot as plt


"""
Problem 3 - Part a
Purpose - Plot trajectories of 16 particles under Lennard-Jone Potential using Verlet algorithm

# Pseudocode
# Define constants - mass, sigma, epsilon, number of particles, total steps, size of time step
# Adapt the code snippet in the instruction
# Define functions
#   acceleration - compute acceleration with the expression derived in Q2.a
#   verlet - compute verlet algorithm 
# Create arrays to hold the positions and velocities of particles
# Save the initial values for positions and velocities
# Run for loop 1000 times (T = 1000)
#   T=0: Compute v(r+0.5h) for the first step
#   T!=0: Apply Verlet algorithm to compute the positions and velocities in the next iteration
# Plot
"""
m = 1  # mass of particle
sigma = 1  # given
epsilon = 1  # given

N = 16 # Number of particles
dt = 0.01 # Time step
T = 1000 # Total steps

Lx = 4.0 # X-limit
Ly = 4.0 # Y-limit

dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)

x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)

xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)

x_initial = xx_grid.flatten() # Create the initial positions of particles
y_initial = yy_grid.flatten()

def acceleration(x, y): # Define a function to calculate acceleration from Lennard-Jones Potential
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

def verlet(dt, x, y, v_x, v_y): # Define a function to compute Verlet algorithm
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

# Arrays to store positions, velocities. Note x_velocity_h and y_velocity_h are to compute kinetic energy in part b
x_position = []
y_position = []
x_velocity = []
y_velocity = []
x_velocity_h = []
y_velocity_h = []

for i in range(N): # Save the initial positions and initial velocities
    x_position.append([x_initial[i]])
    y_position.append([y_initial[i]])
    x_velocity.append([0])
    y_velocity.append([0])
    x_velocity_h.append([0])
    y_velocity_h.append([0])

for t in range(T): # 1000 time steps
    if t == 0: # First step calculates v(t+0.5h)
        for i in range(N):
            particle_x = x_position[i][-1] # Store the latest position of i-th particle
            particle_y = y_position[i][-1]
            x_vel = 0
            y_vel = 0
            for j in range(N): # Loop through the rest of 15 particles
                if i != j: # Exclude the i-th particle
                    particle_j_x = x_position[j][-1]
                    particle_j_y = y_position[j][-1]
                    x_relative = particle_x - particle_j_x
                    y_relative = particle_y - particle_j_y
                    x_acceleration = acceleration(x_relative, y_relative)[0] # Compute Acceleration
                    y_acceleration = acceleration(x_relative, y_relative)[1]
                    x_vel += 0.5 * dt * x_acceleration # Compute v(t+0.5h)
                    y_vel += 0.5 * dt * y_acceleration

            x_velocity[i].append(x_vel) # Save the velocity to the array
            y_velocity[i].append(y_vel)
    else:

        x_position_dt = x_position # Make a copy of the arrays because the original arrays will be modified during the loop
        y_position_dt = y_position
        x_velocity_dt = x_velocity
        y_velocity_dt = y_velocity
        for i in range(N):
            particle_x = x_position_dt[i][-1]   # Store the latest position of i-th particle
            particle_y = y_position_dt[i][-1]
            x_vel = x_velocity_dt[i][-1] # Store the latest velocity of i-th particle
            y_vel = y_velocity_dt[i][-1]
            particle_x_new = particle_x + dt * x_vel # Compute r(t+h)
            particle_y_new = particle_y + dt * y_vel
            x_acceleration = 0
            y_acceleration = 0

            for j in range(N):
                if i != j:  # Exclude the i-th particle
                    particle_j_x = x_position_dt[j][-1]
                    particle_j_y = y_position_dt[j][-1]
                    x_relative = particle_x_new - particle_j_x # The relative distance between i-th and j-th particles
                    y_relative = particle_y_new - particle_j_y

                    x_acceleration += acceleration(x_relative, y_relative)[0] # Compute the acceleration
                    y_acceleration += acceleration(x_relative, y_relative)[1]

            x_position[i].append(particle_x_new) # Store the new position for i-th particle
            y_position[i].append(particle_y_new)
            k_x = dt * x_acceleration # Compute k = h*f(r(t+h), t+h)
            k_y = dt * y_acceleration
            x_velocity[i].append(x_vel + k_x) # Compute v(t + 1.5h)
            y_velocity[i].append(y_vel + k_y)
            x_velocity_h[i].append(x_vel + 0.5*k_x) # Compute v(t+h) for calculation of kinetic energies
            y_velocity_h[i].append(y_vel + 0.5*k_y)


plt.figure(1)
plt.plot(x_position[0], y_position[0], '.', label="Particle 1")
plt.plot(x_position[1], y_position[1], '.', label="Particle 2")
plt.plot(x_position[2], y_position[2], '.', label="Particle 3")
plt.plot(x_position[3], y_position[3], '.', label="Particle 4")
plt.plot(x_position[4], y_position[4], '.', label="Particle 5")
plt.plot(x_position[5], y_position[5], '.', label="Particle 6")
plt.plot(x_position[6], y_position[6], '.', label="Particle 7")
plt.plot(x_position[7], y_position[7], '.', label="Particle 8")
plt.plot(x_position[8], y_position[8], '.', label="Particle 9")
plt.plot(x_position[9], y_position[9], '.', label="Particle 10")
plt.plot(x_position[10], y_position[10], '.', label="Particle 11")
plt.plot(x_position[11], y_position[11], '.', label="Particle 12")
plt.plot(x_position[12], y_position[12], '.', label="Particle 13")
plt.plot(x_position[13], y_position[13], '.', label="Particle 14")
plt.plot(x_position[14], y_position[14], '.', label="Particle 15")
plt.plot(x_position[15], y_position[15], '.', label="Particle 16")
plt.xlabel("x (arbitrary unit)", fontsize=18)
plt.ylabel("y (arbitrary unit)", fontsize=18)
plt.title("Trajectories of 16 particles under Lennard-Jones Potential with T = 1000, dt=0.01", fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left", fontsize=12)


"""
Problem 3 - Part b
Purpose - Confirm the conservation of energy for the system in part a

# Pseudocode
# Create arrays to hold the potential energy and kinetic energy
# Define functions
#   potential: Compute potential energy (Lennard-Jones Potential)
#   kinetic: Compute kinetic energy (0.5mv^2)
# Run for loop 1000 times (T = 1000)
#   Compute potential and kinetic energies using the defined functions
# Plot
"""

potential_energy = []
kinetic_energy = []

def potential(x, y):
    '''
    This is  used to calculate the Lennard-Jones potential

    :param x: relative distance between two particles in x-direction
    :param y: relative distance between two particles in x-direction

    :return: computed value of Lennard-Jones potential
    '''
    r = (x**2 + y**2)**0.5
    return 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

def kinetic(v_x, v_y):
    '''
    This is  used to calculate the kinetic energy 0.5mv^2

    :param v_x: x-component of velocity of a particle
    :param v_y: y-component of velocity of a particle

    :return: computed value of kinetic energy
    '''
    return 0.5 * m * (v_x**2 + v_y**2)

# Calculate Potential and Kinetic energy

for t in range(T): # 1000 time steps
    potential_sum = 0
    kinetic_sum = 0

    for i in range(N):
        particle_x = x_position[i][t] # Store the latest position of i-th particle
        particle_y = y_position[i][t]

        for j in range(N):
            if i != j:
                particle_j_x = x_position[j][t]
                particle_j_y = y_position[j][t]
                x_relative = particle_x - particle_j_x
                y_relative = particle_y - particle_j_y

                potential_sum += 0.5*potential(x_relative, y_relative) # The potential per particle is half Lennard-Jones Potential


        kinetic_sum += kinetic(x_velocity_h[i][t], y_velocity_h[i][t])



    potential_energy.append(potential_sum)
    kinetic_energy.append(kinetic_sum)


plt.figure(2)
plt.plot(np.arange(T), potential_energy,'.', label="Potential Energy")
plt.plot(np.arange(T), kinetic_energy, '.', label="Kinetic Energy")
plt.xlabel("Time Step", fontsize=16)
plt.ylabel("Energy [Arbitrary Unit]", fontsize=16)
plt.title("Potential and Kinetic energy at each time step", fontsize=20)
plt.legend(fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show() # Display a plot

print("Potential Energy")
print(potential_energy)
print("=====================================================")
print("Kinetic Energy")
print(kinetic_energy)



