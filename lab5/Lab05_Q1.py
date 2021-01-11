import numpy as np
import matplotlib.pyplot as plt
import math

"""
Exercise 7.2 Part a
Purpose - Read the data from sunspots.txt and plot a graph of sunspots as a function of time

# Pseudocode
# Place the textfile sunspots.txt into the same folder
# Load the data 
# Store time and the number of sunspots into each array separately
# Plot using pyplot
"""
sunspots_data = np.loadtxt("sunspots.txt")   # Load the data

time = sunspots_data[:,0]    # Store the first column
sunspots = sunspots_data[:,1]    # Store the second column

plt.figure(1)       # Plot of time vs the number of sunspots
plt.plot(time,sunspots,label="The number of sunspots")
plt.xlabel("Time", fontsize=20)
plt.ylabel("The number of sunspots", fontsize=20)
plt.title("The number of sunspots for each month since January 1749", fontsize=20)
plt.legend(fontsize=20)



"""
Exercise 7.2 Part b
Purpose - Plot the power spectrum of the sunspot signal

# Pseudocode
#
"""

def dft(y):  # Custom code of discrete fourier transform
    N = len(y)
    c = np.zeros(N//2+1, complex)
    for k in range(N//2+1):
        for n in range(N):
            c[k] += y[n]*math.exp(-2j*math.pi*k*n/N)
    return c



c_sunspots = np.fft.rfft(sunspots) # Calculate fourier coefficients
k = np.arange(0, len(c_sunspots))

plt.figure(2)   # Power Spectrum of sunspots signal
plt.plot(k, abs(c_sunspots)**2, label="The magnitude squared of the Fourier coefficients")
plt.xlabel("k", fontsize=20)
plt.ylabel("|c$_{k}|^{2}$", fontsize=20)
plt.title("Power Spectrum Of Sunspots Signal", fontsize=20)
plt.legend(fontsize=20)


"""
Exercise 7.4 
Purpose - Perform inverse fourier transform on Dow Jones Industrial Average

# Pseudocode
# Load dow.txt
# Perform Fast Fourier Transform to calculate fourier coefficients
# Set all but the first 10% of the fourier coefficients to zero
# Perform inverse fourier transform and plot
# Set all but the first 2% of the fourier coefficients to zero
# Perform inverse fourier transform and plot
"""

dow_data = np.loadtxt("dow.txt")    # Load dow.txt


c_dow = np.fft.rfft(dow_data)       # Perform Discrete Fourier Transform
set_zeros = np.zeros(len(c_dow))    # Set all but the first ten percent of values to zero

for i in range(int(len(c_dow)*0.1)):
    set_zeros[i] = 1




c_dow_ten = c_dow * set_zeros
c_dow_ten_inverse = np.fft.irfft(c_dow_ten) # Perform inverse fourier transform


set_zeros = np.zeros(len(c_dow))    # Set all but the first two percent of values to zero
for i in range(int(len(c_dow)*0.02)):
    set_zeros[i] = 1

c_dow_two = c_dow * set_zeros
c_dow_two_inverse = np.fft.irfft(c_dow_two) # Perform inverse fourier transform

plt.figure(3)   # Plot (Original + first ten percent + first two percent)
plt.plot(dow_data, label="Original Data")
plt.plot(c_dow_ten_inverse, label="Fast Fourier Transform - Ten percent")
plt.plot(c_dow_two_inverse, label="Fast Fourier Transform - Two percent")
plt.title("Closing Value of the Dow Jones Industrial Average from late 2006 until the end of 2010", fontsize=20)
plt.xlabel("Time", fontsize=20)
plt.ylabel("Dow Jones Industrial Average", fontsize=20)
plt.legend(fontsize=20)
#plt.close()


"""
Exercise 7.6
Purpose - Investigate an additional artifact observed in fourier coefficient graph for dow jones industrial average
          from 2004 until 2008

# Pseudocode

Part a
# Load dow2.txt
# Perform fast fourier transform to calculate fourier coefficients
# Set all but the first 2 percent of the fourier coefficients to zero
# Perform inverse fourier transform and plot

Part b
# Perform discrete cosine transformation
# Set all but the first 2 percent of the fourier coefficients to zero
# Perform inverse fourier transform and plot
"""

dow_data2 = np.loadtxt('dow2.txt') # Load the data

c_dow2 = np.fft.rfft(dow_data2)

set_zeros = np.zeros(len(c_dow2))
for i in range(int(len(c_dow2)*0.02)): # Set all but the first two percent of values to zero
    set_zeros[i] = 1

c_dow2_two = c_dow2 * set_zeros
c_dow2_two_inverse = np.fft.irfft(c_dow2_two) # Perform inverse fourier transform

# Discrete cosine transformation
def dct(y):
    N = len(y)
    y2 = np.empty(2*N,float)
    y2[:N] = y[:]
    y2[N:] = y[::-1]

    c = np.fft.rfft(y2)
    phi = np.exp(-1j*np.pi*np.arange(N)/(2*N))
    return np.real(phi*c[:N])

# Inverse discrete cosine transformation
def idct(a):
    N = len(a)
    c = np.empty(N+1,complex)

    phi = np.exp(1j*np.pi*np.arange(N)/(2*N))
    c[:N] = phi*a
    c[N] = 0.0
    return np.fft.irfft(c)[:N]

c_dow2_cosine = dct(dow_data2)  # Discrete cosine transform
set_zeros = np.zeros(len(c_dow2_cosine))
for i in range(int(len(c_dow2_cosine)*0.02)): # Set all but the first two percent of values to zero
    set_zeros[i] = 1
c_dow2_cosine_two = c_dow2_cosine * set_zeros
c_dow2_cosine_two_inverse = idct(c_dow2_cosine_two) # Inverse discrete cosine transform

plt.figure(4)
plt.plot(dow_data2, label='Original Data')
plt.plot(c_dow2_two_inverse, label="Fast Fourier Transform - Two percent")
plt.plot(c_dow2_cosine_two_inverse, label="Discrete Cosine Transformation - Two percent")
plt.title("Closing Value of the Dow Jones Industrial Average from 2004 until 2008", fontsize=20)
plt.xlabel("Time", fontsize=20)
plt.ylabel("Dow Jones Industrial Average", fontsize=20)
plt.legend(fontsize=20)




"""
Exercise 7.3
Purpose - Find the musical note played by both piano and trumpet when the recording was made

# Pseudocode
# Load the data
# Perform discrete fourier transform on the loaded data to calculate fourier coefficients
# Plot the original data
# Plot the magnitude of the first 10,000 fourier coefficients
"""

piano_data = np.loadtxt('piano.txt')    # Load the data from piano.txt
trumpet_data = np.loadtxt("trumpet.txt")    # Load the data from trumpet.txt

c_piano = np.fft.rfft(piano_data)
c_trumpet = np.fft.rfft(trumpet_data)


plt.figure(5) # Plot the original data of Piano.txt
plt.plot(piano_data, label="Original Data - Piano")
plt.xlabel("# of Sample", fontsize=20)
plt.ylabel("Magnitude", fontsize=20)
plt.title("Waveform of musical notes played on a piano - 100,000 Samples", fontsize=20)
plt.legend(fontsize=20)


plt.figure(6) # Plot the magnitude of the first 10,000 coefficients for Piano
plt.plot(abs(c_piano[:10000]), label="First 10,000 coefficients")
plt.xlabel("# of Sample", fontsize=20)
plt.ylabel("Magnitude of Fourier Coefficients", fontsize=20)
plt.title("Magnitudes of the first 10,000 Fourier coefficients - Piano", fontsize=20)
plt.legend(fontsize=20)


plt.figure(7) # Plot the original data of Trumpet.txt
plt.plot(trumpet_data, label="Original Data - Trumpet")
plt.xlabel("# of Sample", fontsize=20)
plt.ylabel("Magnitude", fontsize=20)
plt.title("Waveform of musical notes played on a trumpet - 23,852 Samples", fontsize=20)
plt.legend(fontsize=20)


plt.figure(8) # Plot the magnitude of the first 10,000 coefficients for Trumpet
plt.plot(abs(c_trumpet[:10000]), label="First 10,000 coefficients")
plt.xlabel("# of Sample", fontsize=20)
plt.ylabel("Magnitude of Fourier Coefficients", fontsize=20)
plt.title("Magnitudes of the first 10,000 Fourier coefficients - Trumpet", fontsize=20)
plt.legend(fontsize=20)

plt.show()