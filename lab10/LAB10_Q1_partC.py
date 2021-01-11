import numpy as np
import matplotlib.pyplot as plt
import random as rd
import matplotlib.animation as animation
import time
"""
Question 1 - Part c
Purpose - Diffusion-limited aggregation until the centre is anchored

# Pseudo code
# Define constants
#   L: Length of lattice
# Create arrays to hold coordinates 
#   x_array: 2D array to hold x-positions (Initial position = (L-1)//2)
#   y_array: 2D array to hold y-positions (Initial position = (L-1)//2)
#   x_anchored: 1D array to hold the x-position of anchored particles
#   y_anchored: 1D array to hold the y-posiiton of anchored particles
#   anchored_board: LxL matrix that stores the location of anchored particles
# Define a function that generates a random number between 0 and 3
#   random.randrange(int)
# Define a function to check hitting the boundary
# Define a function to check hitting the anchored particles
# Use while loop that iterates until the centre of the grid is anchored
#   Call the function to check hitting the boundary
#   Call the function to check hitting the anchored particles
#   Store the position to x_array and y_array
# Define an animation
#  FuncAnimation: frames = length of each list inside x_array and y_array (per particle)
# Run the animation
# Plot the final DLA distribution with x_anchored and y_anchored arrays
"""

start_time = time.time()


L = 151                         # Length of grid
#P = L * L                       # Number of particles
#N = 100000

x_coord = (L-1)//2              # Initial x-position
y_coord = (L-1)//2              # Initial y-position
x_array = []                    # Record all moves in x-direction
y_array = []                    # Record all moves in y-direction
x_anchored = []                 # Record the final position of each particle
y_anchored = []                 # Record the final position of each particle
anchored_board = np.zeros((L, L))   # Board that records the occupied locations

x_array.append([x_coord])       # Set the initial position of a particle
y_array.append([y_coord])

# Function to generate a random number between 0, 1, 2, 3
def random_generator():             # 0 - up, 1 - down, 2 - left, 3 - right
    return rd.randrange(4)

# Check if the particle hit the boundary
# Return True if it did
# Return False otherwise
def check_boundary(x, y):
    if x == 0 or x == L-1 or y == 0 or y == L-1:
        return True

    return False

# Check if the neighboring points are already occupied (anchored_board == 1)
# If occupied, return True.
# Otherwise, return False
def check_anchored(x, y):

    if x == 0:                          # Leftmost
        if anchored_board[x+1][y] == 1:
            return True
        elif anchored_board[x][y-1] == 1:
            return True
        elif anchored_board[x][y+1] == 1:
            return True
    elif x == L-1:                      # Rightmost
        if anchored_board[x - 1][y] == 1:
            return True
        elif anchored_board[x][y - 1] == 1:
            return True
        elif anchored_board[x][y + 1] == 1:
            return True
    elif y == 0:                        # Top
        if anchored_board[x - 1][y] == 1:
            return True
        elif anchored_board[x + 1][y] == 1:
            return True
        elif anchored_board[x][y + 1] == 1:
            return True
    elif y == L-1:                      # Bottom
        if anchored_board[x - 1][y] == 1:
            return True
        elif anchored_board[x + 1][y] == 1:
            return True
        elif anchored_board[x][y - 1] == 1:
            return True
    else:                               # Inner
        if anchored_board[x - 1][y] == 1:
            return True
        elif anchored_board[x + 1][y] == 1:
            return True
        elif anchored_board[x][y - 1] == 1:
            return True
        elif anchored_board[x][y + 1] == 1:
            return True
    return False

num_particle = 0
print("Loop begin")

# Run as long as the center is not occupied
while anchored_board[(L-1)//2][(L-1)//2] != 1:

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

    stick_boundary = check_boundary(new_x, new_y)               # Check if the particle hit the boundary
    stick_anchored = check_anchored(new_x, new_y)               # Check if the position after move has an occupied point in a neighborhood
    original_anchored = check_anchored(x_coord, y_coord)        # Check if the position before move has an occupied point in a neighborhood

    if original_anchored:                                       # If original_anchored == True, record the current position and increase the index of particle
        x_array[num_particle].append(x_coord)
        y_array[num_particle].append(y_coord)
        x_anchored.append(x_coord)
        y_anchored.append(y_coord)
        anchored_board[x_coord][y_coord] = 1                    # Set the current position of the anchored_board to 1
        x_coord = (L - 1) // 2                                  # Initialize the position of the next particle to be the center
        y_coord = (L - 1) // 2

        num_particle += 1                                       # Increase the index of the particle
        x_array.append([x_coord])
        y_array.append([y_coord])

    elif stick_anchored:                                        # If stick_anchored == True, record the next position and increase the index of particle
        x_array[num_particle].append(new_x)
        y_array[num_particle].append(new_y)
        x_anchored.append(new_x)
        y_anchored.append(new_y)
        anchored_board[new_x][new_y] = 1                        # Set the current position of the anchored_board to 1
        x_coord = (L - 1) // 2                                  # Initialize the position of the next particle to be the center
        y_coord = (L - 1) // 2

        num_particle += 1                                       # Increase the index of the particle
        x_array.append([x_coord])
        y_array.append([y_coord])

    elif stick_boundary:                                        # If stick_boundary == True, record the next position and increase the index of particle
        x_array[num_particle].append(new_x)
        y_array[num_particle].append(new_y)
        x_anchored.append(new_x)
        y_anchored.append(new_y)
        anchored_board[new_x][new_y] = 1                        # Set the current position of the anchored_board to 1
        x_coord = (L - 1) // 2                                  # Initialize the position of the next particle to be the center
        y_coord = (L - 1) // 2

        num_particle += 1
        x_array.append([x_coord])
        y_array.append([y_coord])

    else:
        x_array[num_particle].append(new_x)
        y_array[num_particle].append(new_y)
        x_coord = new_x
        y_coord = new_y

print("The number of particles:", len(x_array))

# Run for the number of particles
# In order to see the animation, Make the following changes
# 1) Uncomment the below for loop
# 2) In plt.errorbar, set index_particle+1 to index_particle
# 3) Uncomment myAnimation
"""
for index_particle in range(len(x_array)):
    fig, ax = plt.subplots()

    ax = plt.axes(xlim=(0, L-1), ylim=(0, L-1), xlabel="x", ylabel="y", title="Diffusion-limited aggregation - Until the center is reached")
    plt.grid()

    redDot, = plt.plot([], [], 'ro')                                # Create the red dot

    def init():                                                     # Initialize the position of the red dot
        redDot.set_data([], [])
        return redDot,

    def animate(i):                                                 # Return the move of the red dot
        redDot.set_data(x_array[index_particle][i], y_array[index_particle][i])
        return redDot,
    plt.errorbar(x_anchored[:index_particle+1], y_anchored[:index_particle+1],  fmt='o', color='red')       # Print the anchored particle in each run
#    myAnimation = animation.FuncAnimation(fig, animate, init_func=init, frames=len(x_array[index_particle]), interval=0.001, blit=True, repeat=False)
    fig.canvas.manager.full_screen_toggle()                             # Full Screen
    plt.show(block=False)
    plt.pause(len(x_array[index_particle])*0.001*(10**-3)+0.5)            # Wait 0.5 seconds until the next move
    plt.close()
"""
print("Loop End")
end_time = time.time()

print("The center reached")
print("Run time is ", end_time - start_time, "s")

# Comment the below part and uncomment the above part in order to see the full trajectory
plt.errorbar(x_anchored, y_anchored, fmt='o', color='red')
plt.xlabel("x", fontsize=24)
plt.ylabel("y", fontsize=24)
plt.title("Diffusion-limited aggregation - Until the center is reached", fontsize=24)
plt.xlim(0, L-1)
plt.ylim(0, L-1)
plt.grid()
plt.show()