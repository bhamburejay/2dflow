import json
import subprocess
import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
import pdb


# This strcture defines the inputs that will passed through a .json file to the
# c++ vischyrdo code.  It is relativiely important and its structure should not
# be changed without deliberation. The values can be modified, e.g.
# data['eta_over_s'] = 20.
input_definition = { 
                   "NDOF": 7,   # Number of fields in a cell on the grid. 
                   "NX": 201,
                   "xmin": -60.,
                   "xmax": 60.,
                   "initial_time": 0.,
                   "final_time": 47,
                   "cfl_max": 0.49,
                   "eta_over_s": 1./(4.0*np.pi), 
                   "iofilename": "vischydro.h5",
                   "steps_per_print": 10,
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
data = FixedDataDictionary(input_definition)

# one shoudl keep these options unless you know what you are doing
options = {
    '-ts_type': 'eimex',
    '-ts_max_snes_failures': '5',
    '-pc_type': 'jacobi',
    '-ksp_type': 'bcgs',
    '-ts_adapt_type': 'none', 
}

# Dump the current input data to a json file (default name is
# 'inputs.json'), which serves as input to the vischydro code.
# Additional arguments can be passed to the code on the command
# line by adding them to the options data structure, see option.
#
# The initialdata array is written to an HDF5 file with filename
# data['iofilename'], creating an array initialdata. This is then
# read by the vischydro code. 
def runcode(xarray, initialdata, data, runcommand='./vischydro', inputs='inputs.json', actually_run=True):

    with open(inputs, 'w') as file:
        json.dump(data, file, indent=4)

    # Create an HDF5 file
    iofilename = data['iofilename']
    with h5py.File(iofilename, 'w') as file:
        # Create a dataset in the file and write the array to it
        file.create_dataset('initialdata', data=initialdata) 
        file.create_dataset('x', data=xarray) 
        file.attrs.update(data)
        file.attrs['dt'] = data['cfl_max']*(xarray[1] - xarray[0])

    opts = []
    for kev, val in options.items():
        if val is None or val == '':
            opts += [kev]
        else:
            opts += [kev, val]

    runcommand_array = runcommand.split()
    command = runcommand_array + ['-inputs', inputs] + opts[:] 
    print("Executing command: ")
    print(" ".join(command))
    if actually_run:
        subprocess.call(command)


# Create a grid of x values based on the data dictionary
def xgrid(data):
    # create an array of linear spaced data with NX points between xmin and xmax
    xmin = data['xmin']
    xmax = data['xmax']
    NX = data['NX']
    xarray = np.linspace(xmin, xmax, NX)
    dx = xarray[1] - xarray[0] 
    return xarray, dx

# Create the initial conditions considered by Pretorious and Pandyas fig 1
def ic1(data, A=0.4, delta=0.1, w=25.):
    xarray, dx = xgrid(data) 
    yarray = A*np.exp(-xarray**2/w) + delta
    initialdata = np.zeros((data['NX'], data['NDOF'])) 
    initialdata[:,2] = yarray
    return xarray, initialdata

# Create the step initial conditions considered by Pretorious and Pandyas.
def ic2(data, eps_l=1.0, eps_r=0.1):
    xarray, dx = xgrid(data) 
    yarray = eps_l * np.piecewise(xarray, [xarray < 0, xarray >= 0], [1, 0]) + eps_r * np.piecewise(xarray, [xarray < 0, xarray >= 0], [0, 1])
    initialdata = np.zeros((data['NX'], data['NDOF'])) 
    initialdata[:,2] = yarray
    return xarray, initialdata


