import numpy as np
import matplotlib.pyplot as plt

"""
Problem 3 - Part c
Purpose - Apply periodic boundary conditions to Q3.a
          1) Make particles that exit the domain re-enter
          2) Compute forces due to particles in 8 additional tiles

# Pseudocode
# Define constants - mass, sigma, epsilon, number of particles, total steps, size of time step
# Adapt the code snippet in the instruction
# Define functions
#   acceleration - compute acceleration with the expression derived in Q2.a
#   verlet - compute verlet algorithm 
# Create arrays to hold the positions and velocities of particles
# Create arrays of positions and velocities for particles in 8 additional tiles
# Save the initial values for positions and velocities
# Run for loop 1000 times (T = 1000)
#   T=0: Compute v(r+0.5h) for the first step
#   T!=0: Apply Verlet algorithm to compute the positions and velocities in the next iteration 
          Include calculations of accelerations due to particles in 8 additional tiles
# Plot
"""

m = 1  # mass of particle
sigma = 1  # given
epsilon = 1  # given

N = 16 # Number of particles
dt = 0.01
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

            x_velocity[i].append(x_vel)
            y_velocity[i].append(y_vel)
    else:

        x_position_dt = x_position # Make a copy of the arrays because the original arrays will be modified during the loop
        y_position_dt = y_position
        x_velocity_dt = x_velocity
        y_velocity_dt = y_velocity

        x_position1_dt = np.array(x_position) - 4 # position of mirrored charges in 1st tile
        y_position1_dt = np.array(y_position) + 4

        x_position2_dt = x_position               # position of mirrored charges in 2nd tile
        y_position2_dt = np.array(y_position) + 4

        x_position3_dt = np.array(x_position) + 4 # position of mirrored charges in 3nd tile
        y_position3_dt = np.array(y_position) + 4

        x_position4_dt = np.array(x_position) - 4 # position of mirrored charges in 4nd tile
        y_position4_dt = y_position

        x_position5_dt = np.array(x_position) + 4 # position of mirrored charges in 5nd tile
        y_position5_dt = y_position

        x_position6_dt = np.array(x_position) - 4 # position of mirrored charges in 6nd tile
        y_position6_dt = np.array(y_position) - 4

        x_position7_dt = x_position               # position of mirrored charges in 7nd tile
        y_position7_dt = np.array(y_position) - 4

        x_position8_dt = np.array(x_position) + 4 # position of mirrored charges in 8nd tile
        y_position8_dt = np.array(y_position) - 4

        for i in range(N):
            particle_x = x_position_dt[i][-1] # Store the latest position of i-th particle
            particle_y = y_position_dt[i][-1]
            x_vel = x_velocity_dt[i][-1] # Store the latest velocity of i-th particle
            y_vel = y_velocity_dt[i][-1]
            particle_x_new = particle_x + dt * x_vel # Compute r(t+h)
            particle_y_new = particle_y + dt * y_vel
            x_acceleration = 0
            y_acceleration = 0

            for j in range(N):
                if i != j: # Exclude the i-th particle
                    particle_j_x = x_position_dt[j][-1]
                    particle_j_y = y_position_dt[j][-1]
                    x_relative = particle_x_new - particle_j_x # The relative distance between i-th and j-th particles
                    y_relative = particle_y_new - particle_j_y

                    x_acceleration += acceleration(x_relative, y_relative)[0] # Compute the acceleration
                    y_acceleration += acceleration(x_relative, y_relative)[1]

                particle_k1_x = x_position1_dt[j][-1]
                particle_k1_y = y_position1_dt[j][-1]
                x_relative_k1 = particle_x_new - particle_k1_x
                y_relative_k1 = particle_y_new - particle_k1_y
                x_acceleration += acceleration(x_relative_k1, y_relative_k1)[0]
                y_acceleration += acceleration(x_relative_k1, y_relative_k1)[1]

                particle_k2_x = x_position2_dt[j][-1]
                particle_k2_y = y_position2_dt[j][-1]
                x_relative_k2 = particle_x_new - particle_k2_x
                y_relative_k2 = particle_y_new - particle_k2_y
                x_acceleration += acceleration(x_relative_k2, y_relative_k2)[0]
                y_acceleration += acceleration(x_relative_k2, y_relative_k2)[1]

                particle_k3_x = x_position3_dt[j][-1]
                particle_k3_y = y_position3_dt[j][-1]
                x_relative_k3 = particle_x_new - particle_k3_x
                y_relative_k3 = particle_y_new - particle_k3_y
                x_acceleration += acceleration(x_relative_k3, y_relative_k3)[0]
                y_acceleration += acceleration(x_relative_k3, y_relative_k3)[1]

                particle_k4_x = x_position4_dt[j][-1]
                particle_k4_y = y_position4_dt[j][-1]
                x_relative_k4 = particle_x_new - particle_k4_x
                y_relative_k4 = particle_y_new - particle_k4_y
                x_acceleration += acceleration(x_relative_k4, y_relative_k4)[0]
                y_acceleration += acceleration(x_relative_k4, y_relative_k4)[1]

                particle_k5_x = x_position5_dt[j][-1]
                particle_k5_y = y_position5_dt[j][-1]
                x_relative_k5 = particle_x_new - particle_k5_x
                y_relative_k5 = particle_y_new - particle_k5_y
                x_acceleration += acceleration(x_relative_k5, y_relative_k5)[0]
                y_acceleration += acceleration(x_relative_k5, y_relative_k5)[1]

                particle_k6_x = x_position6_dt[j][-1]
                particle_k6_y = y_position6_dt[j][-1]
                x_relative_k6 = particle_x_new - particle_k6_x
                y_relative_k6 = particle_y_new - particle_k6_y
                x_acceleration += acceleration(x_relative_k6, y_relative_k6)[0]
                y_acceleration += acceleration(x_relative_k6, y_relative_k6)[1]

                particle_k7_x = x_position7_dt[j][-1]
                particle_k7_y = y_position7_dt[j][-1]
                x_relative_k7 = particle_x_new - particle_k7_x
                y_relative_k7 = particle_y_new - particle_k7_y
                x_acceleration += acceleration(x_relative_k7, y_relative_k7)[0]
                y_acceleration += acceleration(x_relative_k7, y_relative_k7)[1]

                particle_k8_x = x_position8_dt[j][-1]
                particle_k8_y = y_position8_dt[j][-1]
                x_relative_k8 = particle_x_new - particle_k8_x
                y_relative_k8 = particle_y_new - particle_k8_y
                x_acceleration += acceleration(x_relative_k8, y_relative_k8)[0]
                y_acceleration += acceleration(x_relative_k8, y_relative_k8)[1]

            x_position[i].append(np.mod(particle_x_new, Lx)) # Store the new position for i-th particle
            y_position[i].append(np.mod(particle_y_new, Ly)) # Check if the particle crossed the wall
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
plt.xlabel("x (arbitrary unit)", fontsize=16)
plt.ylabel("y (arbitrary unit)",fontsize=16)
plt.title("Trajectory of particle 1 under Lennard-Jones Potential with T = 1000, dt=0.01 (Periodic Boundary Condition)", fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left", fontsize=12)
plt.show()
