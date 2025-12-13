import json
from collections import UserDict
import subprocess
import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
import pdb


# inputs for the vischydro code
_input_data = {
    "Vischydro": {
        "nx": 80,
        "ny": 90,
        "xmin": -8.0,
        "xmax": 8.0,
        "ymin": -9.0,
        "ymax": 9.0,

        # These are defaults and should not normally be changed
        "ndof": 9,
        "cfl_max": 0.8,
        "is_bjorken": True,
        "highest_order_term_only": False
    },
    "VischydroMain": {
        "print_frequency": 10,
        "t_start": 0.0,
        "t_end": 47.0,
        "dt_max": 1.0,
        "run_name": "DFTest",
        "initial_conditions_filename": "DFTest_ic.h5",
        "eta_by_s": 1./(4.*np.pi),
        "zeta_by_s": 0.0
    }
}


class FixedData(UserDict):
    def __init__(self, initialdata):
        # super().__init__(initialdata)
        self.data = initialdata

    def __setitem__(self, key, value):
        """ Navigate the nested dictionary using the provided keystring
        checking that each key exists before setting the value to v. The format
        of keystring is 'key1/key2/key3' for nested dictionaries. """
        keystring = key.split('/')
        name = keystring[-1]
        d = self.data
        for i, key in enumerate(keystring[:-1]):
            d = d[key]
        if d.get(name) is None:
            raise KeyError("Key '{}' not found in dictionary.".format(name))
        d[name] = value

    def __getitem__(self, key):
        """ Navigate the nested dictionary using the provided keystring
        checking that each key exists before getting the value to v. The format
        of keystring is 'key1/key2/key3' for nested dictionaries. """
        keystring = key.split('/')
        name = keystring[-1]
        d = self.data
        for i, key in enumerate(keystring[:-1]):
            d = d[key]
        if d.get(name) is None:
            print(keystring)
            raise KeyError("Key '{}' not found in dictionary.".format(name))
        return d[name]


input_data = FixedData(_input_data)

# Create a Vischydro of xy values based on the input_data dictionary


def xygrid(input_data):
    # create an array of linear spaced input_data with NX points between xmin and xmax
    xmin = input_data['Vischydro']['xmin']
    xmax = input_data['Vischydro']['xmax']
    ymin = input_data['Vischydro']['ymin']
    ymax = input_data['Vischydro']['ymax']
    NX = input_data['Vischydro']['nx']
    xarray = np.linspace(xmin, xmax, NX, endpoint=False)
    NY = input_data['Vischydro']['ny']
    yarray = np.linspace(ymin, ymax, NY, endpoint=False)
    dx = xarray[1] - xarray[0]
    dy = yarray[1] - yarray[0]
    # create a meshgrid
    xarray, yarray = np.meshgrid(xarray, yarray, indexing='xy')
    print("xygrid: ", "Created Vischydro coordinate grid based on the input data with nx =", NX, " ny =", NY)
    return xarray, yarray, dx, dy


def initialdata_grid(input_data):
    # Creates a grid to hold the intiial data
    nx = input_data['Vischydro']['nx']
    ny = input_data['Vischydro']['ny']
    ndof = input_data['Vischydro']['ndof']
    initialdata = np.zeros(
        (input_data['Vischydro']['ny'], input_data['Vischydro']['nx'], input_data['Vischydro']['ndof']))
    print("initialdata_grid:", "Created Vischydro coordinate grid based on the input data with nx =",
          nx, " ny =", ny, " ndof =", ndof)
    return initialdata


# one shoudl keep these options unless you know what you are doing
options = {}

# potential options
#    '-ts_type': 'eimex',
#    '-ts_max_snes_failures': '5',
#    '-pc_type': 'jacobi',
#    '-ksp_type': 'bcgs',
#    '-ts_adapt_type': 'none',

def runcode(initialdata, input_data, runcommand='./VischydroMain.exe', actually_run=True, petsc_args=''):
    """ Run VicschydroMain.exe

    Dump the current input input_data to a json file (default name is
    'RUNNAME.json'), which serves as input to the vischydro code.
    Additional arguments can be passed to the code on the command
    line by adding them to the options input_data structure, see option.
    
    The initialdata array is written to an HDF5 file with filename
    input_data['icfilename'], creating an array initialdata. This is then
    read by the vischydro code. """

    # Create an HDF5 file
    icfilename = input_data['VischydroMain/run_name'] + '_ic.h5'

    input_data['VischydroMain/initial_conditions_filename'] = icfilename
    with h5py.File(icfilename, 'w') as file:
        # Create a dataset in the file and write the array to it
        file.create_dataset('initialdata', data=initialdata)

    inputs = input_data['VischydroMain/run_name'] + '.json'
    with open(inputs, 'w') as file:
        json.dump(input_data.data, file, indent=4)

    opts = []
    for kev, val in options.items():
        if val is None or val == '':
            opts += [kev]
        else:
            opts += [kev, val]
    if petsc_args:
        opts += petsc_args.split()

    runcommand_array = runcommand.split()
    command = runcommand_array + ['-input', inputs] + opts[:]
    print("Executing command: ")
    print(" ".join(command))
    if actually_run:
        subprocess.call(command)


# Swap the Vischydro dimensions for x and y
def swap_xy(input_data):
    # swap nx and ny
    nx = input_data['Vischydro']['nx']
    ny = input_data['Vischydro']['ny']
    input_data['Vischydro/nx'] = ny
    input_data['Vischydro/ny'] = nx
    # swap xmin and ymin
    xmin = input_data['Vischydro']['xmin']
    ymin = input_data['Vischydro']['ymin']
    input_data['Vischydro/xmin'] = ymin
    input_data['Vischydro/ymin'] = xmin
    # swap xmax and ymax
    xmax = input_data['Vischydro']['xmax']
    ymax = input_data['Vischydro']['ymax']
    input_data['Vischydro/xmax'] = ymax
    input_data['Vischydro/ymax'] = xmax
