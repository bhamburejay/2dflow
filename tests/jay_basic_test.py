import vischydro
import numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse
import pprint
import sys

runcommand  = 'mpiexec -n 4 ../build/tests/VischydroMain.exe'
   
#runcommand = 'mpiexec -n 4 ./VischydroMain.exe'

def getnames2d(fourpietabys, gaussian_const, gaussian_amplitude, gaussian_width):
    return 'simple_bj2d_example'

def runbj2d(fourpietabys, gaussian_const, gaussian_amplitude, gaussian_width, tmin, tmax, petsc_args='', actually_run=True):
    # These are the inputs that must be passed to the vischydro code through an inputfile in the json format.  You can modify these as needed to change the grid size, tstart and stop times, etc.
    print("Creating example input_data for Vischydro")

    A = gaussian_amplitude 
    delta = gaussian_const
    w = gaussian_width
    fourpietabys = fourpietabys

    name = getnames2d(fourpietabys, delta, A, w)

    # Here will run an example in cartesian
    vischydro.input_data["Vischydro/nx"] = 60 
    vischydro.input_data["Vischydro/xmin"] = -15.0
    vischydro.input_data["Vischydro/xmax"] = 15.0
    vischydro.input_data["Vischydro/ny"] = 60
    vischydro.input_data["Vischydro/ymin"] = -15.0
    vischydro.input_data["Vischydro/ymax"] = 15.0


    vischydro.input_data["VischydroMain/eta_by_s"] = fourpietabys/(4.0*np.pi)
    vischydro.input_data["VischydroMain/zeta_by_s"] =0.0
    vischydro.input_data["VischydroMain/t_start"] = tmin
    vischydro.input_data["VischydroMain/t_end"] = tmax
    vischydro.input_data["VischydroMain/run_name"] = name
    vischydro.input_data["VischydroMain/print_frequency"] = 5
    
    pprint.pprint(vischydro.input_data)

    # Given the initial conditions we can create some initial data
    initialdata, dx, dy = ic1(vischydro.input_data, A=A, delta=delta, sigma=w)

    # pass command line options to vischydro
    
    if actually_run:
        vischydro.runcode(initialdata, vischydro.input_data, runcommand=runcommand, actually_run = actually_run, petsc_args=petsc_args)
    # plot_contours(vischydro.input_data, name)
    plot_tauT00_vs_x_at_fixed_y(vischydro.input_data, name, y_target=2.0)
    plt.savefig(f"{name}_tauT00_slice.png")
    print(f"Plot saved to {name}_tauT00_slice.png")

def plot_tauT00_vs_x_at_fixed_y(input_data, name, y_target=2.0):
    try:
        times = np.loadtxt('{}_grid_t.txt'.format(name))
        # Ensure times is iterable if single line
        if times.ndim == 1 and times.shape[0] == 2: # [step, time]
             times = times.reshape(1, 2)
    except OSError:
        print("Could not load times file.")
        return
        
    print("Times loaded from file: ", times)

    n_steps = times.shape[0]
    # Pick up to 5 times evenly spaced
    num_plots = min(5, n_steps)
    indices = np.linspace(0, n_steps-1, num_plots, dtype=int)
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    with h5py.File(name + '_grid.h5', 'r') as file:
        # Coordinates shape of dataset "coordinates" is usually (ny, nx, 2) or similar
        # Let's verify shape. Usually it's stored once or per step? 
        # In Vischydro.cpp: VecView(run->coordinates, H5viewer); 
        # HDF5 structure: 'coordinates' dataset.
        coords = file['coordinates'][:]
        # coords shape: (ny, nx, 2) where 2 is (x, y)
        
        # Find index for y closest to y_target
        # y is the second coordinate, so coords[:, 0, 1] gives y-values column
        y_col = coords[:, 0, 1] 
        iy = (np.abs(y_col - y_target)).argmin()
        actual_y = y_col[iy]
        print(f"Taking slice at index {iy}, y = {actual_y:.3f}")
        
        x_slice = coords[iy, :, 0]
        
        for idx in indices:
            time_val = times[idx, 1] 
            
            # E is index 0 in solution (time, y, x, dof)
            E_slice = file['solution'][idx, iy, :, 0]

            # Calculate tau * T00 => time_val * E_slice
            y_vals = time_val * E_slice
            
            ax.plot(x_slice, y_vals, label=f'$\\tau={time_val:.2f}$')

    ax.set_xlabel('x [fm]')
    ax.set_ylabel(r'$\tau T^{00}$')
    ax.legend()
    ax.set_title(r'$\tau T^{{00}}$ vs x at y $\approx$ {:.1f}'.format(actual_y))
    ax.grid(True)
    
# Create the initial conditions considered by Pretorious and Pandyas fig 1
def ic1(input_data, A=0.4, delta=0.1, sigma=4.):
    # Get the xy grid 
    xarray, yarray, dx, dy = vischydro.xygrid(input_data)

    zarray = A*np.exp(-xarray**2/(2.*sigma**2) - yarray**2/(2.*(sigma*1.4)**2)) + delta

    # creates a grid with the right dimensions for the vischydro code
    # the initial values are all zero 
    initialdata = vischydro.initialdata_grid(input_data) 

    # set the energy density initial condition slots (3rd index = 3) The third
    # index is the energy density. 
    initialdata[:,:, 3] = zarray 

    return initialdata, dx, dy

def plot_contours(input_data, name):
    times  = np.loadtxt('{}_grid_t.txt'.format(name))
    print("Times loaded from file: ", times)

    to_plot = [0, times.shape[0]//3, 2*times.shape[0]//3, -1]
    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.flatten()
    with h5py.File(name + '_grid.h5', 'r') as file:
        x = file['coordinates'][:, :,0]
        y = file['coordinates'][:, :,1]

        for i, it in enumerate(to_plot):
            solution = file['solution'][it,  :, :, 0] # energy density is the 4th variable
            etabys = input_data['VischydroMain/eta_by_s']
            axs[i].set_title('t = {:.2f} $4\\pi\\eta/s$ = {}'.format(times[it,0],4.0*np.pi*etabys))
            axs[i].contourf(x, y, solution, levels=50)
            

def example_bj2d(petsc_args='', run=False):
    fourpietabys = [2.]
    gaussian_const = 0.000001

    Ce = 2. * (3.0**2 -1.) * np.pi**2 /30.0 
    tau0T0 = 1.
    tau0 = 1.0
    T0 = tau0T0 / tau0
    e0 = Ce * T0**4
    R = 3.5 

    Kn = max(1./(4 * np.pi) * (4.0 * fourpietabys[0]) / (3.0 * tau0 * T0), 1.e-10)
    print("Initial inverse Knudsen number for 4pi eta/s = {} is 1/Kn = {}".format(fourpietabys[0], 1/Kn))
    
    gaussian_amplitude = e0
    gaussian_width = R
    tmax = 12.0

    runbj2d(fourpietabys[0], gaussian_const, gaussian_amplitude, gaussian_width, tau0, tmax, petsc_args=petsc_args, actually_run=run)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run example bj2d viscous hydrodynamics simulation and plot results for a simple case.')
    parser.add_argument('--dont_run', help='Do not run the example_bj2d test, just plot the results from a previous run.', action='store_true')
    args, unknown = parser.parse_known_args()

    example_bj2d(petsc_args='', run=not args.dont_run)

        
