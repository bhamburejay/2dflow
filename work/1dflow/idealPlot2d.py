import h5py
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Load data from HDF5 file
with h5py.File('idealOutput2d.h5', 'r') as f:
    initial = f['initialdatain'][:]
    final = f['finaldata'][:]

# Assume shape is (NY, NX, NDOF)
NY, NX, NDOF = initial.shape
x = np.linspace(-60, 60, NX)
y = np.linspace(-60, 60, NY)

# Energy density index (match 1D code: e is index 2)
e_idx = 2
e_init = initial[:, :, e_idx]
e_final = final[:, :, e_idx]

# 1D slice along x at center y
center_y = NY // 2
plt.figure(figsize=(8, 4))
plt.plot(x, e_init[center_y, :], label='Initial $e$ (y=0)')
plt.plot(x, e_final[center_y, :], label='Final $e$ (y=0)')
plt.xlabel('x')
plt.ylabel('Energy density $e$')
plt.title('Energy Density Slice at y=0')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# 1D slice along y at center x
center_x = NX // 2
plt.figure(figsize=(8, 4))
plt.plot(y, e_init[:, center_x], label='Initial $e$ (x=0)')
plt.plot(y, e_final[:, center_x], label='Final $e$ (x=0)')
plt.xlabel('y')
plt.ylabel('Energy density $e$')
plt.title('Energy Density Slice at x=0')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
plt.grid(True)
