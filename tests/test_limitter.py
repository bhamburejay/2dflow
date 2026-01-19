import numpy as np

# This is a script to test the numerical derivative calculation done by test_limitter.cpp in the DFHydro code. First we generate the test data ./test_limitter and then we read in the data, from test_limitter_slope.dat,  and compare to the analytical derivative.

data = np.loadtxt('test_limitter_slope.dat')
x, f, df = data[:,0], data[:,1], data[:,2]
dx = x[1] - x[0]
sigma = 1.0
df_analytical = dx * (-x / sigma**2 * np.exp(-x**2 / (2*sigma**2)))

# The first and last points have bad numerical derivatives due to boundary effects. 
error = np.abs(df[1:-1] - df_analytical[1:-1])
print("Max error:", np.max(error))
print("RMS error:", np.sqrt(np.mean(error**2)))

import matplotlib.pyplot as plt

# Plot the function f(x)
plt.figure(figsize=(10, 6))
plt.plot(x, f, label='f(x)', color='blue')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Function f(x)')
plt.legend()
plt.grid()
plt.show()

# Plot the numerical and analytical derivatives
plt.figure(figsize=(10, 6))
plt.plot(x[1:-1], df[1:-1], '.', label='Numerical Derivative', color='red', linestyle='--')
plt.plot(x[1:-1], df_analytical[1:-1], label='Analytical Derivative', color='green')
plt.xlabel('x')
plt.ylabel('df/dx')
plt.title('Derivatives of f(x)')
plt.legend()
plt.grid()
plt.show()

# Plot the error between numerical and analytical derivatives
plt.figure(figsize=(10, 6))
plt.plot(x[1:-1], error, label='Error', color='purple')
plt.xlabel('x')
plt.ylabel('Error')
plt.title('Error between Numerical and Analytical Derivatives')
plt.legend()
plt.grid()
plt.show()
