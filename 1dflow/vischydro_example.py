import numpy as np
import matplotlib.pyplot as plt
import h5py
import runvischydro as rv

# The code will run with the default options. The options can be changed by updating the options dictionary. But you could just leave the defaults
rv.options.update({'-ts_exact_final_time': 'INTERPOLATE'})
print("The options passed to vischydro are:")
print(rv.options)



# Set the values of other paraemters (see runvischydro.py)
# The data associate array in the rv module contains
# all the parameters which are passed to the program. 
# Here we override the default value of eta over s
rv.data['eta_over_s'] = 3./(4.0*np.pi) 

# This data file is the one used to communicate the initial conditions and final conditions
rv.data['iofilename'] = 'vischydro_fig1.h5'

# print out the values of the parameters passed to the code
print("The input passed to vischydro are:")
print(rv.data)

# construct the initial data in an array. Read the code.
# this returns the xcoordinates and the initial energy density
xarray, edensity = rv.ic1(rv.data, A=0.48, delta=0.12, w=25 )
#print(xarray, edensity)

# Run the code, with the initial data and the chosen options.
# The inputs are read from the file myinputs.json
#
# On Derek's machine runcommand is mpiexec-mpich-clang17 -n 1 ./vischyro'
rv.runcode(xarray, edensity, rv.data, runcommand='mpiexec -n 1 ./vischydro -ksp_monitor_error -log_view', inputs='myinputs.json')

# Open the HDF5 file which is used for input and output
with h5py.File('vischydro_fig1.h5', 'r') as file:
    # Read the 'finaldata' dataset
    finaldata = file['finaldata'][:]
    initialdatain = file['initialdatain'][:]

    plt.plot(xarray, finaldata[:, 2], xarray, finaldata[:, 0], '--')
    ax = plt.gca()
    ax.set_xlim(-60, 60)
    ax.set_ylim(0.09, 0.25)
    plt.legend([r'$\epsilon$', r'$T^{tt}$'])
    plt.xlabel('x')
    plt.title(r'$\eta$/s = {}/4$\pi$'.format(rv.data['eta_over_s']*4.0*np.pi)) 
    plt.show()
        

