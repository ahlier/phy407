import numpy as np

'''
This program use monte carlo method to calculate the volume of n-dimensinal 
hypersphere

pseudo-code:
- write a function that check if a given point is inside the hypersphere, if yes
return True, else return False
- create an array of random number ranging from -1 to 1, size dimension X sampling 
number, and put check how many point are inside hypersphere
- calculate volume of hypersphere and error using equation given
'''

# set constant
N = 1000000  # number of sampling points
dim = 10  # dimension


def f(x):
    # this function takes in array of sampling points, return array of
    # True/False, True if point is inside hypersphere, False if not
    r_square = np.zeros(len(x[0]))

    for i in range(len(x)):
        for j in range(len(x[i])):
            r_square[j] += x[i][j]**2
    return r_square < 1


x = (np.random.random((dim, N))-0.5)*2
# random will only generate value of 0-1, this will make random values range
# from -1 to 1

fx = f(x)
I = 2 ** dim / N * sum(fx)  # volume of hypersphere

var = sum(fx ** 2) / N - (sum(fx) / N) ** 2
sigma = 2 ** dim * np.sqrt(var / N)  # error
print('I = {} +/- {}'.format(I, sigma))

