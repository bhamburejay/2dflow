import numpy as np

data = np.loadtxt('test_slope.dat')
x, f, df = data[:,0], data[:,1], data[:,2]
sigma = 1.0
df_analytical = -x / sigma**2 * np.exp(-x**2 / (2*sigma**2))
error = np.abs(df - df_analytical)
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
plt.plot(x, df, label='Numerical Derivative', color='red', linestyle='--')
plt.plot(x, df_analytical, label='Analytical Derivative', color='green')
plt.xlabel('x')
plt.ylabel('df/dx')
plt.title('Derivatives of f(x)')
plt.legend()
plt.grid()
plt.show()

# Plot the error
plt.figure(figsize=(10, 6))
plt.plot(x, error, label='Error', color='purple')
plt.xlabel('x')
plt.ylabel('Error')
plt.title('Error between Numerical and Analytical Derivatives')
plt.legend()
plt.grid()
plt.show()