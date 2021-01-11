import numpy as np
import matplotlib.pyplot as plt

'''
This program will turn a blur picture into a clear one

pseudo code:
- read data from blur.txt file
- create point spread function, and a matrix using point spread function
- use fourier transformation on matrices in blur.txt and point spread function
- divide blur picture with point spread function, and check if the point spread 
function is less than epsilon, if, yes, don't change the value. 
'''
sigma = 25  # given, used in point spread function

epsilon = 10 ** -4
# if the fourier of spread function is less than epsilon, we don't change the
# fourier of blurred picture

file = np.loadtxt('blur.txt', float)  # read the picture
y_position, x_position = file.shape  # get the size of photo

# Show the blur picture
plt.imshow(file)
plt.title('Blurred Image')
plt.show()


def point_spread(x, y):
    # The point spread function is given, x, y are x,y position
    return np.e**(- ( x ** 2 + y ** 2 ) / (2 * sigma ** 2))


# update the value of each point using point spread function
point_spread_array = np.ones([ x_position, y_position ], float)
for i in range(x_position):
    for j in range(y_position):
        # The point spread function should have x,y value centered around the
        # origin, the following will count from origin to the right/up, then
        # from origin to the left/down
        point_spread_array[i, j] = point_spread( (j + x_position / 2) % x_position - x_position / 2, \
                                                 (i + y_position / 2) % y_position - y_position / 2)

# use fourier transformation on blur picture and point spread function
# fft2 will normalize the function, so we don't need to calculate the coefficient
fourier_blur = np.fft.fft2(file)
fourier_point_spread = np.fft.fft2(point_spread_array)


# The following divide blur picture fourier transformation with spread function
# fourier transformation
unblurred_fourier = np.empty([ x_position, y_position // 2 + 1], complex)
for i in range(y_position // 2 + 1):
    for j in range(x_position):
        # check if the fourier transformation of spread function is larger than
        # epsilon, keep the value of point spread function if fourier of point
        # spread function is less than epsilon
        if abs(fourier_point_spread[j, i]) < epsilon:
            unblurred_fourier[j, i] = fourier_blur[j, i]
        else:
            unblurred_fourier[j, i] = fourier_blur[j, i] / (fourier_point_spread[j, i])

# get the inverse fourier transformation of unblurred_fourier
plt.imshow(np.fft.irfft2(unblurred_fourier))
plt.title('Clear Image')
plt.gray()
plt.show()


# Q2_b
# create a matrix using point spread function
spread_graph = np.ones((x_position, y_position))
for i in range(y_position):
    for j in range(x_position):
        spread_graph[i, j] = point_spread( (j + x_position / 2) % x_position - x_position / 2, \
                                                 (i + y_position / 2) % y_position - y_position / 2)

# multiple by the picture matrix to get the density matrix
spread_graph *= file
plt.imshow(spread_graph)
plt.title('Density plot')
plt.show()

