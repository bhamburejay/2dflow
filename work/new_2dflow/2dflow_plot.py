import numpy as np
import matplotlib.pyplot as plt
import h5py

def plot_initial(tag='2dflow_test'):
    # Load the initial state
    with h5py.File(tag + '_initial.h5', 'r') as f:
        coordinates = f['coordinates'][:]
        # 'output' contains [nx, ny, nvar], where nvar=number of hydro variables
        data2 = f['output'][:, :, 0]  # Energy density is usually the first variable

        X = coordinates[:, :, 0]
        Y = coordinates[:, :, 1]

        plt.figure()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.contourf(X, Y, data2, levels=100)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Initial Energy Density')
        plt.colorbar(label='Energy Density')

def plot_central_evolution(tag='2dflow_test'):
    # Load the time-evolving solution
    with h5py.File(tag + '_grid.h5', 'r') as f:
        coordinates = f['coordinates'][:]
        solution = f['solution'][:]  # shape: [nt, nx, ny, nvar]
        # Find the center indices
        ix = coordinates.shape[0] // 2
        iy = coordinates.shape[1] // 2
        # Extract energy density at the center for all times
        e = solution[:, ix, iy, 0]  # [nt]
    # Load time grid from ascii file
    tgrid = np.loadtxt(tag + '_grid_t.txt')
    t = tgrid[:, 0]  # First column is time

    plt.figure()
    plt.plot(t, e, '.', label='Numerical')
    plt.xlabel('t')
    plt.ylabel('Energy Density at Center')
    plt.title('Central Energy Density vs Time')
    plt.legend()
    plt.grid(True)

# Run the plotting functions
plot_initial()
plot_central_evolution()
plt.show()
