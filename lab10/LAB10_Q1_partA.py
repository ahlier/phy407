import numpy as np
import matplotlib.pyplot as plt
import random as rd
import matplotlib.animation as animation

"""
Question 1 - Part a
Purpose - Demonstrate the Brownian motion with a random walk (Exercise 10.3)

# Pseudo code
# Define constants
#   L: Length of lattice
#   N: Number of steps
# Create arrays to hold coordinates 
#   x_array: 1D array to hold x-positions (Initial position = (L-1)//2)
#   y_array: 1D array to hold y-positions (Initial position = (L-1)//2)
# Define a function that generates a random number between 0 and 3
#   random.randrange(int)
# Define a function to check hitting the boundary
# Use for loop that iterates N number of moves
#   Call the function to check hitting the boundary
#   Use while loop to keep making move until the particle doesn't go over boundary
#   Store the position to x_array and y_array
# Define an animation
#   FuncAnimation: N number of frames
# Run the animation
"""

L = 101                         # Length of Lattice
N = 5000                        # Number of steps

x_coord = (L-1)//2              # Initial x-position
y_coord = (L-1)//2              # Initial y-position
x_array = [x_coord]             # Record all moves in x-direction
y_array = [y_coord]             # Record all moves in y-direction

# Function to generate a random number between 0, 1, 2, 3
def random_generator():             # 0 - up, 1 - down, 2 - left, 3 - right
    return rd.randrange(4)

# Check if the particle hit the boundary
# Return True if it did
# Return False otherwise
def check_boundary(x, y):
    if 0 <= x <= 100 and 0 <= y <= 100:
        return False
    else:
        return True

# Run for N number of steps
for i in range(N):
    outside_boundary = True
    while(outside_boundary):
        move = random_generator()
        new_x = x_coord
        new_y = y_coord
        if move == 0:   # up
            new_y += 1
        elif move == 1: # down
            new_y -= 1
        elif move == 2: # left
            new_x -= 1
        elif move == 3: # right
            new_x += 1

        outside_boundary = check_boundary(new_x, new_y)
    x_coord = new_x                 # Set the x-coordinate to the next move
    y_coord = new_y                 # Set the xy-coordinate to the next move
    x_array.append(x_coord)         # Set the next move to x_array
    y_array.append(y_coord)         # Set the next move to y_array


fig, ax = plt.subplots()

#ax = plt.axis([0, L-1, 0, L-1])                     # Axis Range is from 0 to 100
ax = plt.axes(xlim=(0, L-1), ylim=(0, L-1), xlabel="x", ylabel="y", title="Brownian motion (1 particle) - 5000 steps")
plt.grid()

redDot, = plt.plot(x_array[0], y_array[0], 'ro')    # Create the red dot
def animate(i):                                     # Return the move of the red dot
    redDot.set_data(x_array[i], y_array[i])
    return redDot,

# Create an animation for N number of steps
myAnimation = animation.FuncAnimation(fig, animate, frames=np.arange(1, N), interval=1, blit=True, repeat=False)
fig.canvas.manager.full_screen_toggle()             # Full Screen
plt.show(block=False)
plt.pause(N * 1 * (10 ** -3) + 1)                   # Wait one seconds after the animation is over
plt.close()
