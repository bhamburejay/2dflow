import vischydro
import numpy as np
import h5py
import matplotlib.pyplot as plt
import pprint

def example_test():
    # These are the inputs that must be passed to the vischydro code through an inputfile in the json format.  You can modify these as needed to change the grid size, tstart and stop times, etc.
    print("Creating example input_data for Vischydro")
    pprint.pprint(vischydro.input_data)

    # Here will run an example in cartesian
    vischydro.input_data["Vischydro"]["is_bjorken"] = False

    # Given the initial conditions we can create some initial data
    initialdata, dx, dy = ic1(vischydro.input_data)

    # (1) Writes the initialdata to an h5 file that can be read by vischydro
    # called input_data['VischydroMain']['run_name'] + '_ic.h5'
    # 
    # (2) This also updates the input_data dictionary to set the initial_conditions_filename key
    # 
    # (3) Then writes the input_data to a jsonfile that can be read by vischydro
    # called input_data['VischydroMain']['run_name'] + '.json'
    #
    # This just writes the command that would be used to run vischydro, it does not actually run it.
    vischydro.runcode(initialdata, vischydro.input_data, runcommand='mpiexec-mpich-clang17 -n 1 ./VischydroMain.exe', actually_run = False)

    # Now run the code  
    vischydro.runcode(initialdata, vischydro.input_data, runcommand='mpiexec-mpich-clang17 -n 1 ./VischydroMain.exe', actually_run = True)
    


# Create the initial conditions considered by Pretorious and Pandyas fig 1
def ic1(input_data, A=0.4, delta=0.1, w=25., swapxy=False):
    # Get the xy grid 
    xarray, yarray, dx, dy = vischydro.xygrid(input_data)
    if swapxy:
        zarray = A*np.exp(-yarray**2/w) + delta
    else:
        zarray = A*np.exp(-xarray**2/w) + delta

    # creates a grid with the right dimensions for the vischydro code
    # the initial values are all zero 
    initialdata = vischydro.initialdata_grid(input_data) 

    # set the energy density initial condition slots (3rd index = 3) The third
    # index is the energy density. Indices 4, 5 are the velocity components, ux
    # and uy. See VischydroNode.h One should only set the indices 3, 4, 5. All
    # others are filled by the code. 
    initialdata[:,:, 3] = zarray 

    return initialdata, dx, dy


def plot_default_output_slice(input_data, filename, itime = 0, swapxy=False):
    with h5py.File(filename, 'r') as file:
        finaldata = file['solution'][itime,:,:,:]
        coordinates = file['coordinates'][:]
        # print (finaldata.shape)
        # print (coordinates.shape)
        if swapxy:
            ny, nx, ndof = finaldata.shape
            #xVischydro = coordinates[ny//2,:,0]
            yVischydro = coordinates[:, nx//2,1]
            energy_density = finaldata[:,nx//2,3]
            plt.plot(yVischydro, energy_density, '.', label='{}'.format(itime))
        else:
            ny, nx, ndof = finaldata.shape
            Vischydro = coordinates[ny//2,:,0]
            #yVischydro = coordinates[:, nx//2,1]
            energy_density = finaldata[ny//2,:,3]
            plt.plot(xVischydro, energy_density, '.', label='{}'.format(itime))


# Plots a central slice throught the input_data file
def plot_default_output(input_data, filename, swapxy=False):
    with h5py.File(filename, 'r') as file:
        finaldata = file['output'][:,:,:]
        coordinates = file['coordinates'][:]
        if swapxy:
            ny, nx, ndof = finaldata.shape
            #xVischydro = coordinates[:,nx//2,0]
            yVischydro = coordinates[:, nx//2,1]
            energy_density = finaldata[:, nx//2, 3]
            plt.plot(yVischydro, energy_density, '.', label='{}'.format(filename))
            return
        else:
            ny, nx, ndof = finaldata.shape
            xVischydro = coordinates[ny//2,:,0]
            #yVischydro = coordinates[:, nx//2,:,1]
            energy_density = finaldata[ny//2,:,3]
            plt.plot(xVischydro, energy_density, '.', label='{}'.format(filename))
            return


if __name__ == "__main__":
    example_test()

    idata = vischydro.input_data
    run_name = vischydro.input_data['VischydroMain']['run_name'] 
    plot_default_output(idata, run_name + '_initial.h5', swapxy=False)
    plot_default_output(idata, run_name + '_final.h5', swapxy=False)
    plt.legend()
    plt.show()
