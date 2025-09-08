import h5py
import numpy as np
import matplotlib.pyplot as plt

# Load data from HDF5 file
with h5py.File('final_idealhydro2d.h5', 'r') as f:
    initial = f['initialdata'][:]
    final = f['finaldata'][:]

# Assume shape is (NY, NX, 9) for 9 dof
NY, NX, _ = initial.shape
x = np.linspace(-60, 60, NX)
y = np.linspace(-60, 60, NY)
X, Y = np.meshgrid(x, y)

# Energy density

e_init = initial[:, :, 3]  # e is at index 3
e_final = final[:, :, 3]

# 2D contour plots
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.contourf(X, Y, e_init, levels=50, cmap='viridis')
plt.colorbar(label='Initial Energy Density $e$')
plt.title('Initial Energy Density (2D)')
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(1, 2, 2)
plt.contourf(X, Y, e_final, levels=50, cmap='viridis')
plt.colorbar(label='Final Energy Density $e$')
plt.title('Final Energy Density (2D)')
plt.xlabel('x')
plt.ylabel('y')

plt.tight_layout()
plt.show()

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

# Velocity field (ux, uy)

ux_init = initial[:, :, 4]
uy_init = initial[:, :, 5]
ux_final = final[:, :, 4]
uy_final = final[:, :, 5]

# Streamplot for initial velocity over initial energy density
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.contourf(X, Y, e_init, levels=50, cmap='viridis', alpha=0.7)
plt.streamplot(X, Y, ux_init, uy_init, color='k', density=1.2, linewidth=1)
plt.title('Initial Velocity Field (Streamplot)')
plt.xlabel('x')
plt.ylabel('y')

# Streamplot for final velocity over final energy density
plt.subplot(1, 2, 2)
plt.contourf(X, Y, e_final, levels=50, cmap='viridis', alpha=0.7)
plt.streamplot(X, Y, ux_final, uy_final, color='k', density=1.2, linewidth=1)
plt.title('Final Velocity Field (Streamplot)')
plt.xlabel('x')
plt.ylabel('y')

plt.tight_layout()
plt.show()
