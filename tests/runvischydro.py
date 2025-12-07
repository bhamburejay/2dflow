import json
import subprocess
import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
import pdb

# inputs for the vischydro code
input_definition = {
    "Vischydro": {
        "nx": 201,
        "ny": 7,
        "ndof": 9, 
        "xmin": -60.0,
        "xmax": 60.0, 
        "ymin": -3.0,
        "ymax": 3.0, 
        "cfl_max": 0.49
    },
    "VischydroMain": {
        "print_frequency" : 10, 
        "t_start": 0.0,
        "t_end": 47.0,
        "run_name" : "DFTest", 
        "initial_conditions_filename": "DFTest_ic.h5"
    }
}

class FixedDataDictionary(dict):
    """
    A dictionary with a fixed set of keys
    """

    def __init__(self, dictionary):
        dict.__init__(self)
        for key in dictionary.keys():
            dict.__setitem__(self, key, dictionary[key])

    def __setitem__(self, key, item):
        if key not in self:
            raise KeyError("The key {} is not defined".format(key))
        dict.__setitem__(self, key, item)

# 
# Creaate a dictionary with fixed keys.
#
input_data = FixedDataDictionary(input_definition)

# Swap the Vischydro dimensions for x and y
def swap_xy_Vischydro(input_data):
    # swap nx and ny
    nx = input_data['Vischydro']['nx']
    ny = input_data['Vischydro']['ny']
    input_data['Vischydro']['nx'] = ny
    input_data['Vischydro']['ny'] = nx
    # swap xmin and ymin
    xmin = input_data['Vischydro']['xmin']
    ymin = input_data['Vischydro']['ymin']
    input_data['Vischydro']['xmin'] = ymin
    input_data['Vischydro']['ymin'] = xmin
    # swap xmax and ymax
    xmax = input_data['Vischydro']['xmax']
    ymax = input_data['Vischydro']['ymax']
    input_data['Vischydro']['xmax'] = ymax
    input_data['Vischydro']['ymax'] = xmax


# one shoudl keep these options unless you know what you are doing
options = {
    '-ts_type': 'ssp',
    '-ts_max_snes_failures': '5',
    '-pc_type': 'jacobi',
    '-ksp_type': 'bcgs',
    '-ts_adapt_type': 'none', 
}

# Dump the current input input_data to a json file (default name is
# 'RUNNAME.json'), which serves as input to the vischydro code.
# Additional arguments can be passed to the code on the command
# line by adding them to the options input_data structure, see option.
#
# The initialdata array is written to an HDF5 file with filename
# input_data['icfilename'], creating an array initialdata. This is then
# read by the vischydro code. 
def runcode(initialdata, input_data, runcommand='./VischydroMain.exe', actually_run=True):

    # Create an HDF5 file
    icfilename = input_data['VischydroMain']['run_name'] + '_ic.h5'
    input_data['VischydroMain']['initial_conditions_filename'] = icfilename

    with h5py.File(icfilename, 'w') as file:
        # Create a dataset in the file and write the array to it
        file.create_dataset('initialdata', data=initialdata)

    inputs = input_data['VischydroMain']['run_name'] + '.json'
    with open(inputs, 'w') as file:
        json.dump(input_data, file, indent=4)

    opts = []
    for kev, val in options.items():
        if val is None or val == '':
            opts += [kev]
        else:
            opts += [kev, val]

    runcommand_array = runcommand.split()
    command = runcommand_array + ['-input', inputs]  + opts[:] 
    print("Executing command: ")
    print(" ".join(command))
    if actually_run:
        subprocess.call(command)


# Create a Vischydro of xy values based on the input_data dictionary
def xyVischydro(input_data):
    # create an array of linear spaced input_data with NX points between xmin and xmax
    xmin = input_data['Vischydro']['xmin']
    xmax = input_data['Vischydro']['xmax']
    ymin = input_data['Vischydro']['ymin']
    ymax = input_data['Vischydro']['ymax'] 
    NX = input_data['Vischydro']['nx']
    xarray = np.linspace(xmin, xmax, NX)
    NY = input_data['Vischydro']['ny']
    yarray = np.linspace(ymin, ymax, NY)
    dx = xarray[1] - xarray[0]
    dy = yarray[1] - yarray[0]
    # create a meshgrid
    xarray, yarray = np.meshgrid(xarray, yarray, indexing='xy')         
    print(xarray.shape, yarray.shape)
    # print(xarray)
    # print(yarray)
    return xarray, yarray, dx, dy

# Create the initial conditions considered by Pretorious and Pandyas fig 1
def ic1(input_data, A=0.4, delta=0.1, w=25., swapxy=False):
    xarray, yarray, dx, dy = xyVischydro(input_data)
    if swapxy:
        zarray = A*np.exp(-yarray**2/w) + delta
    else:
        zarray = A*np.exp(-xarray**2/w) + delta

    initialdata = np.zeros((input_data['Vischydro']['ny'], input_data['Vischydro']['nx'], input_data['Vischydro']['ndof'])) 

    # set the energy density initial condition slots (3rd index = 3)
    initialdata[:,:, 3] = zarray

    # make a contour plot of the initial conditions
    # plt.contourf(xarray, yarray, initialdata[:,:, 3], levels=100)
    # plt.colorbar()
    # plt.title("Initial Conditions")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.savefig("initial_conditions.png")
    # plt.show()
    # print("xarray shape:", xarray.shape)
    # print("dx =", dx, xarray[0, 1] - xarray[0, 0], " dy =", dy, yarray[1, 0] - yarray[0, 0])
    
    return initialdata, dx, dy

def plot_default_output_slice(input_data, filename, itime = 0, swapxy=False):
    with h5py.File(filename, 'r') as file:
        finaldata = file['solution'][itime,:,:,:]
        coordinates = file['coordinates'][:]
        print (finaldata.shape)
        print (coordinates.shape)
        ny, nx, ndof = finaldata.shape
        xVischydro = coordinates[ny//2,:,0]
        yVischydro = coordinates[ny//2,:,1]
        energy_density = finaldata[ny//2,:,3]
        plt.plot(xVischydro, energy_density, '.', label='{}'.format(itime))


# Plots a central slice throught the input_data file
def plot_default_output(input_data, filename, swapxy=False):
    with h5py.File(filename, 'r') as file:
        finaldata = file['output'][:,:,:]
        coordinates = file['coordinates'][:]
        print (finaldata.shape)
        print (coordinates.shape)
        ny, nx, ndof = finaldata.shape

        if swapxy:
            xVischydro = coordinates[:,nx//2,1]
            yVischydro = coordinates[:,nx//2,0]
            energy_density = finaldata[:,nx//2,3]
            plt.plot(xVischydro, energy_density, '.', label='{}'.format(filename))
            return
        else:
            xVischydro = coordinates[ny//2,:,0]
            yVischydro = coordinates[ny//2,:,1]
            energy_density = finaldata[ny//2,:,3]
            plt.plot(xVischydro, energy_density, '.', label='{}'.format(filename))
            return

if __name__ == "__main__":
    swapxy= False
    if swapxy:
        swap_xy_Vischydro(input_data)
    initialdata, dx, dy = ic1(input_data, swapxy=swapxy)
    runcode(initialdata, input_data, runcommand='mpiexec-mpich-clang17 -n 1 ./VischydroMain.exe', actually_run=True)

    run_name = input_data['VischydroMain']['run_name'] 
    plot_default_output(input_data, run_name + '_initial.h5', swapxy=swapxy)
    plot_default_output(input_data, run_name + '_final.h5', swapxy=swapxy)
    plt.legend()
    plt.show()
