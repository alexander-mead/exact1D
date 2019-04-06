# Add my own libraries to the path
#import sys
#sys.path.append('/Users/Mead/Physics/library/python')

# Import statements
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Initial white space
print()

# Input file
infile = 'data/output.dat'
print('Input file: ', infile)

# Read the data
data = np.loadtxt(infile)

# Write data info to screen
print('Data size:', np.size(data))
print('Data shape:', np.shape(data))
Nsheet = (np.size(data[0,:])-1)//2
Ntime = np.size(data[:,0])
print('Number of sheets:', Nsheet)
print('Number of time steps:', Ntime)
print()

# Make a big array of all the position coordinates
#pos = np.zeros(Nsheet*Ntime)
pos = []
for i in range(Nsheet):
    pos.extend(data[:,2*i+1])
pos = np.abs(pos) # So that we are considering 'r' rather than 'x'

print('Position array size:', np.size(pos))
print('Position array shape:', np.shape(pos))
print()

# Model profile
rs = 0.1 # Scale radius
a = 10.  # Amplitude
p = 1.   # Power law
def f(r,rs,a,p):
    return a/(1.+(r/rs)**p)

# Create histogram bins
xmin = 1e-3
xmax = 1e0
nx = 128
bins = np.logspace(np.log10(xmin),np.log10(xmax),nx) # Histogram bin edges
x = np.zeros(np.size(bins)-1) # Histogram bin centres
for i in range(np.size(x)):
    x[i]=0.5*(bins[i]+bins[i+1])

# Fit model to data
rho, _ = np.histogram(pos,bins,density=True)
popt, _ = curve_fit(f, xdata=x[0:100], ydata=rho[0:100], p0=[rs,a,p])
rs=popt[0]
a=popt[1]
p=popt[2]
print('Fitting rs:', rs)
print('Fitting a:', a)
print('Fitting p:', p)
print()

# Make the plot
plt.hist(pos,bins,density=True)
plt.gca().set_xscale("log")
plt.gca().set_yscale("log")
plt.plot(x,f(x,rs,a,p),label='power-law',color='black')
#plt.plot(x,rho)
plt.xlabel(r'$r$')
plt.ylabel(r'$\rho(r)$')
plt.ylim((1e-1,2e1))
plt.draw()
plt.pause(0.01)
input("<Hit Enter To Close>")
print()
