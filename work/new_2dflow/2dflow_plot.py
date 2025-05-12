import numpy as np 
import matplotlib.pyplot as plt
import h5py

def plot_initial(tag='2dflow_test'):
    # Load the data in initial.h5 the dataset is called output
    with h5py.File(tag + '_initial.h5', 'r') as f:
        coordinates = f['coordinates'][:]
        data2 = f['output'][:,:,0] 
        # construct a meshgrid from coordinates which
        # will be used to plot the data
        X = coordinates[:,:,0]
        Y = coordinates[:,:,1]
        plt.figure()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.contourf(X, Y, data2)

def plot_central_evolution(tag='2dflow_test'):
    # Load the data in grid.h5 the dataset is called output
    with h5py.File(tag + '_grid.h5', 'r') as f:
        coordinates = f['coordinates'][:]
        print(coordinates.shape)
        ix = coordinates.shape[0]//2 
        iy = coordinates.shape[1]//2
        e = f['solution'][:,ix,iy, 0] 
        tgrid = np.loadtxt(tag + '_grid_t.txt')
        t = tgrid[:,0]
        plt.figure()
        plt.gca().set_xlabel('t')
        plt.gca().set_ylabel('E')
        plt.plot(t, e, '.')
        # plt.plot(t, e[0]*(t[0]/t)**(4./3.))  # Theoretical decay (disabled)

plot_initial()
plot_central_evolution()
plt.show()
    
