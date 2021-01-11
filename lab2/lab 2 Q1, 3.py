import numpy as np
import matplotlib.pyplot as plt

'''
For this program, we need to show that for python, forward difference 
approximation works best when the width of step is 1*10^-8.

pseudo code:
- firstly, write out the equations needed for this program
- then run a for-loop, increase the step width from 10^-16 to 10^0, increase the 
width by 10 for each iteration. Also, Calculate and record the error for each step. 
- Lastly, make a plot of log(step) and log(error), and find the minimum value. 
'''
x = 0.5  # x position of the graph

def f(x):
    # this is the function given
    return np.e**(-x**2)

def forward_diff_f(x, h):
    # this is the derivative of function f at point x using forward different
    # method. This function takes in x position and step width, and output the
    # slope of f at x.
    return (f(x+h)-f(x))/h


def first_derivative_f(x):
    # This is the first derivative of function f, takes in x position and output
    # the derivative at x
    return -2*x*np.e**(-x**2)


def second_derivative_f(x):
    # This is the second derivative of function f, takes in x position and
    # output the second derivative at x
    return -2*np.e**(-2*x**2)+4*x**2*np.e**(-x**2)


file = open('error.txt', 'w')  # store value for step and error
file.write('step\terror\n')
step = []  # used to stored the width of steps
abs_err_forward = []
# used to stored the error of using forward difference scheme
for i in range(17):
    # Each iteration will increase the step width 10 times. And also it will
    # record the step width and error of each iteration
    h = 10**(-16)*10**i
    step.append(h)
    abs_err_forward.append(abs(first_derivative_f(x)-forward_diff_f(x, h)))
    file.write('{}\t{}\n'.format(h, abs_err_forward[-1]))

file.close()
# plot the graph of log(step width) vs log(error)
plt.plot(np.log(step), np.log(abs_err_forward), label='forward difference')
plt.title('Change in Log of Error Magnitude with Respect to Change in Step \n'
          'Width Using Forward Difference Scheme')
plt.xlabel('log of step width')
plt.ylabel('log of error magnitude')
plt.legend()
plt.show()


index = abs_err_forward.index(min(abs_err_forward))
# Finding the step that has minimum error
print('Forward difference scheme has a minimum error when the step size is '
      '{}'.format(step[index]))

# Q1_d  -----------------------------------------------------------
'''
This program is similar to previous one, but replace the forward difference
scheme with central difference scheme
'''
def central_difference_f(x, h):
    # it finds the derivative of function f at point x using central difference
    # scheme. This function takes in x position and step width, and output the
    # slope of f at x.
    return (f(x+h/2)-f(x-h/2))/h


abs_err_central = []  # stores the error when using central difference scheme
for i in range(17):
    # Each iteration will increase the step width 10 times. And also it will
    # record the error of each iteration
    h = 10**(-16)*10**i
    abs_err_central.append(abs(first_derivative_f(x)-central_difference_f(x, h)))

# make a graph of log(step) vs log(error)
plt.plot(np.log(step), np.log(abs_err_forward), label='forward difference')
plt.plot(np.log(step), np.log(abs_err_central), label='central difference')
plt.title('Change in Log of Error Magnitude with Respect to Change in Step \n'
          'Width Using Central Difference Scheme')
plt.xlabel('log of step width')
plt.ylabel('log of error magnitude')
plt.legend()
plt.show()

index = abs_err_central.index(min(abs_err_central))
# Finding the step that has minimum error
print('Forward difference scheme has a minimum error when the step size is '
      '{}'.format(step[index]))

# Q3 -------------------------------------------------------------
# Q3_b
f = 1  # focal length, in meter
slit_separation = 20*10**(-6)  # unit in meter
num_slit = 10  # number of slit
grating_width = slit_separation*num_slit  # width of the grating, in meter
wavelength = 500*10**(-9)  # unit meter
screen_width = 0.1  # unit in meter
N = 1000  # number of step
h = grating_width/N  # width of each step, h is around 10^-8,
# which has the least error from result of Q1
x = np.linspace(-0.05, 0.05, 500)  # observed points on screen
alpha = np.pi/slit_separation


def q(u):
    # this is intensity transmission function, it takes in distance from the
    # center on the grating
    return (np.sin(alpha*u))**2

# Q3_c -------------------------------------------------------------
'''
This program will produce the light intensity graph on a screen based on the 
diffraction grating equation entered.

pseudo code:
- first write out the diffraction grating equation
- then divide the screen into many points, and run a for-loop through all these
points, for each iteration, do an integral of the diffraction grating equation 
and store the values.
- Lastly, with the values we got from integral, plot a graph a position vs 
intensity. Also, make a density graph 
'''

def integrand(x, u):
    # This function can be used to calculate the intensity of light with given
    # grating width u and position on the screen x
    return np.sqrt(q(u)) * np.e**(2j * np.pi * x * u / (wavelength * f))


intensity = [] # stored the intensity of light at every observed point on screen
for j in x:
    # trapezoidal rule is used to calculated the derivative, each iteration
    # will add the light intensity at a point on screen.
    central_difference = 0.5*h*(integrand(j, -grating_width/2)+
                                integrand(j, grating_width/2))
    for i in range(1, N-1):
        central_difference += integrand(j, i*h-grating_width/2)*h
    intensity.append(abs(central_difference)**2)


plt.plot(x, intensity, 'o', label='Data point')
plt.xlabel('position (m)')
plt.ylabel('intensity')
plt.title('Light Intensity on the Screen Using Equation \n'
          'q(u)=sin(alpha*u)^2')
plt.legend()
plt.show()


# Q3_d  -------------------------------------------------------

intensity_array = np.zeros([100, 500])  # create a density plot of size 100 by 500
for i in range(100):
    # add the intensity to each column of the graph
    intensity_array[i, :] = intensity
# make the graph
plt.imshow(intensity_array)
plt.gray()
clb = plt.colorbar()
clb.ax.set_title('intensity', pad = 15)
plt.title('Density Plot of Diffraction Grating Using Equation \n'
          'q(u)=sin(alpha*u)^2')
plt.show()


# Q3_e part 1 --------------------------------------------------
def q_2(u):
    # modify function q, adding new term to it
    return (np.sin(alpha*u))**2*(np.sin(alpha/2)*u)**2


def integrand_2(x, u):
    # same as integrand function, but uses q_2 instead
    return np.sqrt(q_2(u))*np.e**(2j*np.pi*x*u / (wavelength * f))


intensity_2 = []  # stores intensity
for j in x:
    # light intensity is calculated using trapezoidal rule, and each
    # iteration will store the light intensity on the screen
    central_difference_2 = 0.5*h*(integrand_2(j, -grating_width/2)+
                                  integrand_2(j, grating_width/2))
    for i in range(1, N-1):
        central_difference_2 += integrand_2(j, i*h-grating_width/2)*h
    intensity_2.append(abs(central_difference_2)**2)

# make the graph
plt.plot(x, intensity_2, 'o', label='Data point')
plt.xlabel('position (m)')
plt.ylabel('intensity')
plt.title('Light Intensity on the Screen Using Equation \n'
          'q(u)=sin(alpha*u)^2*sin(beta*u)^2')
plt.legend()
plt.show()


intensity_array = np.zeros([100, 500])  # density graph will have size 100x500
for i in range(100):
    # each iteration will add intensity value to column
    intensity_array[i, :] = intensity_2
# make the graph
plt.imshow(intensity_array)
plt.gray()
clb = plt.colorbar()
clb.ax.set_title('intensity', pad=15)
plt.title('Density Plot of Diffraction Grating Using Equation \n'
          'q(u)=sin(alpha*u)^2*sin(beta*u)^2')
plt.show()


# Q3_e part 2
grating_width = 10**(-6)*(20+10+60)  # in meter
N = 90  # number of step, it's taken to be 90, so that h is around 10^-8
h = grating_width/N  # width of each step


intensity_3 = []
for j in x:
    # since the area between slits does not transmit light, so the intensity
    # contribution is 0 in those area, so only calculate the area of slits
    central_difference_3 = 0.5*h*(integrand_2(j, -grating_width/2)+integrand_2(j, grating_width/2))
    for i in range(0, 20):
        central_difference_3 += integrand_2(j, i*h-grating_width/2)*h

    for i in range(80, 90):
        central_difference_3 += integrand_2(j, i*h-grating_width/2)*h

    intensity_3.append(abs(central_difference_3) ** 2)


# make the graph
plt.plot(x, intensity_3, 'o', label='Data point')
plt.xlabel('position (m)')
plt.ylabel('intensity')
plt.title('Light Intensity on the Screen with Two \n'
          'Slits of Different Size')
plt.legend()
plt.show()


intensity_array = np.zeros([100, 500])  # density graph will have size 100x500
for i in range(100):
    # each iteration will add intensity value to column
    intensity_array[i, :] = intensity_3
# make the graph
plt.imshow(intensity_array)
plt.gray()
clb = plt.colorbar()
clb.ax.set_title('intensity', pad=15)
plt.title('Density Plot of Diffraction Grating with 2 Slits of Different Size')
plt.show()
