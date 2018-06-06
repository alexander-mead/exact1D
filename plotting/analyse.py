import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from itertools import cycle
import matplotlib.animation as animation

#Flags
do_trajectory = False
do_virial = False
do_phase = False
do_phase_animation = True
do_sheets = False
do_histogram_xt = False
do_histogram_xv = False

#Read in the positional and velocity data
print('Loading n-body data')
data=np.loadtxt('output.dat')
print('Done')
print()

#Check that the input data looks sensible
print('Data shape:', data.shape)
print('Data type:', type(data))
print()
n=(len(data[0,:])-1)//2
N=len(data[:,0])
print('Number of sheets:', n)
print('Number of time steps:', N)
print()

#Isolate the time steps
t=data[:,0]
print('Length of time array:', len(t))
print('Type of time array:', type(t))
print('Shape of time array:', t.shape)
print('Time array:', t)
tmin=t[0]
tmax=t[N-1]
print('Starting time:', tmin)
print('Final time:', tmax)
print()

#Isolate the position and velocity data
x=data[:,1:2*n:2]
v=data[:,2:2*n+1:2]

#Check that position and velocity look sensible
print('Length of position array:', len(x))
print('Type of position array:', type(x))
print('Shape of position array:', x.shape)
print('Length of velocity array:', len(v))
print('Type of velocity array:', type(v))
print('Shape of velocity array:', v.shape)
print()

#Read in the energy data
print('Loading energy data')
energy=np.loadtxt('energy.dat')
print('Done')
print()

#Check the energy data is sensible and read into arrays with nice names
t2=energy[:,0]
t_equal=np.array_equal(t,t2)
print('Are time arrays equal? ', np.array_equal(t,t2))
if(t_equal == False): exit()
Kinetic=energy[:,1]
Potential=energy[:,2]
Energy=energy[:,3]
Virial=energy[:,4]

if(do_trajectory):

    #Plot the trajectory of all sheets: x(t)
    plt.plot(x,t)
    #plt.plot(t,v)
    plt.xlabel('Position')
    plt.ylabel('Time')
    plt.xlim((-1,1))
    plt.ylim((tmin,tmax))
    plt.draw()
    plt.pause(0.01)
    input('')
    plt.close()

if(do_virial):

    #Plot the virial condtion as a function of time: K/V (t)
    plt.axhline(0.5,c='black')
    plt.plot(t,Virial)
    plt.xlabel('Time')
    plt.ylabel('Kinetic / Potential energies')    
    plt.xlim((tmin,tmax))
    plt.ylim((0,2))
    plt.draw()
    plt.pause(0.01)
    input('')
    plt.close()

if(do_phase):

    #Plot the phase trajectories v(x) for all sheets
    #plt.plot(x[:,11],v[:,11])
    plt.plot(x,v)
    plt.xlabel('x')
    plt.ylabel('v')
    plt.draw()
    plt.pause(0.01)
    input('')
    plt.close()

if(do_phase_animation):

    # Create new Figure and an Axes which fills it.
    fig = plt.figure(figsize = (7, 7))

    cycol = cycle('bgrcmk')

    #scat = ax.scatter(x[0,:],v[0,:])
    plt.xlim((-1,1))
    plt.xlabel('x')
    plt.ylim((-2,2))
    plt.ylabel('v')
    scat = plt.scatter(x[0,:],v[0,:],s=10,c=np.random.random(n))  

    #Function to update plot
    def update(i):
        #print(i, x[i,0], v[i,0])
        data = np.column_stack((x[i,:],v[i,:])) #Column_stack takes in tuple of 'vectors' to stack
        scat.set_offsets(data) #Data must be a 2D array here
        return scat,

    animation = animation.FuncAnimation(fig, update, frames=N, interval=10, repeat=False)
    plt.show()
    plt.pause(0.01)
    input('')
    plt.close()
    exit()

if(do_sheets):

    #Make the figure (axes, canvass etc.)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.xlim(-1,1)
    plt.ylim(0,1)
    ax.set_yticklabels([])
    i=0
    #plt.plot(data[i,:],np.zeros_like(data[i,:]),'x')
    for j in range(n):
        plt.axvline(data[i,j])

    #Routine to update the plot (does not need to return anything)
    def update_plot(i):
        x=L*np.sin(vs[i,1])
        y=-L*np.cos(vs[i,1])
        data=np.stack((x,y),axis=-1)
        im.set_offsets(data)
        line.set_xdata((0.,x))
        line.set_ydata((0.,y))

    #Generate the animation
    #movie=animation.ArtistAnimation(fig,frames,interval=10,repeat=False)

    plt.draw()
    plt.pause(0.01)
    input('')
    plt.close()

if(do_histogram_xt):

    xflat=x.ravel()
    print('xflat:', xflat)
    print('type(xflat):', type(xflat))
    print('len(xflat):', len(xflat))
    print()

    tflat=np.ravel(np.tile(t,(N-1,1)),order='F')
    print('flat:', tflat)
    print('type(tflat):', type(tflat))
    print('len(tflat):', len(tflat))
    print()

    nx=100
    plt.hist2d(xflat,tflat,bins=[nx,N],cmap='gray',normed=True,cmax=0.1)
    plt.xlabel('Position')
    plt.ylabel('Time')
    plt.show()
    plt.pause(0.01)
    input('')
    plt.close()
    exit()

if(do_histogram_xv):

    def frames_append(array,ax):
        frames.append((ax.imshow(array,interpolation='nearest'),))

    #Make the figure (axes, canvass etc.)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    #Create an empty list of all image frames
    frames = []

    print('1')
    phase=np.histogram2d(x[0,:],v[0,:],bins=(100,100))#,range=((-1.,1.),(-2.,2.)))
    print('2')
    frames_append(phase,ax)
    print('3')

    #Generate the animation
    movie=animation.ArtistAnimation(fig,frames,interval=50,repeat=False)
    plt.draw()
    plt.pause(0.01)
    input("<Hit any key to close>")
    plt.close()
