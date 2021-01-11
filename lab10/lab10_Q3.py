import numpy as np
import matplotlib.pyplot as plt

'''
This program uses two methods; importance sampling and mean value method to 
approximate the integral of two equations.

pseudo-code:
- write out the integrating function and probability distribution function
- write a function for mean value method, this function will break the integrating
length into N points, and take the middle value of each two points, multiply by 
the width of these two points and sum them all up.
- write a function for importance sampling, this function will randomly generate
N points. Subbing in these points into integrating function to get weighted 
average of the function. Then sum them up.
- repeat the calculation for mean value method and importance sampling 100 times,
and store the results
- make a histogram of integral values for these two methods. 
'''


N = 10000  # number of points

# integrating range
starting = 0
ending = 1


def f(x):
    # integrating function
    return x**-0.5/(1+np.exp(x))

def p(x):
    # probability distribution
    return 1/(2*np.sqrt(x))


def importance_sampling(N):
    # This function will approximate the integral of function f(x) using N
    # sampling points
    x = np.random.random(N)**2

    def g(x):
        return f(x)/p(x)

    I = sum(g(x))/ N

    return I


def mean_value(N):
    # This function will approximate the integral of function f(x) using N points
    x = np.arange(starting, ending, (ending-starting)/N)
    y = f(x+0.5*(ending-starting)/N)
    return sum(y)*(ending-starting)/N


Is = []  # stored the integral value using importance sampling
Mv = []  # stored the integral value using mean value method

for i in range(100):
    # do the calculation 100 times
    Is.append(importance_sampling(N))
    Mv.append(mean_value(N))

# take the mean of 100 calculations
I_importance = np.mean(Is)
I_mean = np.mean(Mv)
print('The integral calculated using mean value method is {}'.format(I_mean))
print('and Importance sampling give {}'.format(I_importance))

# plot the frequency of integral value
plt.hist(Is, 10, range=[0.8, 0.88], label='Importance Sampling')
plt.xlabel('Integral value')
plt.ylabel('Frequency')
plt.title('Frequency of Integral Value Between 0.8 and 0.88 in 100 Runs \n'
          'Using Importance Sampling Method')
plt.legend()
plt.show()

plt.hist(Mv, 10, range=[0.8, 0.88], label='Mean Value')
plt.xlabel('Integral value')
plt.ylabel('Frequency')
plt.title('Frequency of Integral Value Between 0.8 and 0.88 in 100 runs\n'
          'Using Mean Value Method')
plt.legend()
plt.show()


# part b-----------------------------------------------------------
# part b is the same as part a, but change the integrating function
def f_b(x):
    return np.exp(-2*abs(x-5))


def p_b(x):
    return np.exp((-(x-5)**2)/2)/(np.sqrt(2*np.pi))


def importance_sampling_b(N):
    # The distribution of sampling points should follow gaussian distribution
    x = np.array([np.random.normal(5, 1) for j in range(N)])

    def g_b(x):
        return f_b(x)/p_b(x)

    I = sum(g_b(x))/ N

    return I


def mean_value_b(N):
    x = np.arange(starting, ending, (ending-starting)/N)
    y = f_b(x+0.5*(ending-starting)/N)
    return sum(y)*(ending-starting)/N



starting = 0
ending = 10

Is = []
Mv = []

for i in range(100):
    Is.append(importance_sampling_b(N))
    Mv.append(mean_value_b(N))

I_importance = np.mean(Is)
I_mean = np.mean(Mv)
print('The integral calculated using mean value method is {}'.format(I_mean))
print('and Importance sampling give {}'.format(I_importance))

plt.hist(Is, 10, range=[0.95, 1.05], label='Importance Sampling')
plt.xlabel('Integral value')
plt.ylabel('Frequency')
plt.title('Frequency of Integral Value Between 0.95 and 1.05 in 100 runs\n'
          'Using Importance Sampling')
plt.legend()
plt.show()

plt.hist(Mv, 10, range=[0.95, 1.05], label='Mean Value')
plt.xlabel('Integral value')
plt.ylabel('Frequency')
plt.title('Frequency of Integral Value Between 0.95 and 1.05 in 100 runs\n'
          'Using Mean Value Method')
plt.legend()
plt.show()
