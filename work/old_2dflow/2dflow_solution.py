import numpy as np
import matplotlib.pyplot as plt


file_path = '/Users/space/Desktop/Xhydro/2dflow/build/2dflow_solution.txt'
# Read data from solution.txt
data = np.loadtxt(file_path)  # Load the data from the file
t = data[:, 0]  # Extract the time values (first column)
E = data[:, 1]  # Extract the solution values (second column)

# Plot E(t) vs t
plt.figure(figsize=(10, 6))
plt.plot(t, E, label="E(t)", color="blue", linewidth=2)
plt.xlabel("Time (t)", fontsize=14)
plt.ylabel("Solution (E)", fontsize=14)
plt.title("Solution of the ODE: E(t) vs t", fontsize=16)
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()


# Show the plot
plt.show()