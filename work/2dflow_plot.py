# Note: This file has a code for both 2d and 3d plots. The 3d plots are still work in progress
# 2d works fine

import h5py
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

################## 2D plot #########################

# Read HDF5 data
with h5py.File('energy_out.h5', 'r') as f:
    dataset = f['Energy'][:]
    energy_data = dataset[:, :, 0]  # First component is energy

# NOTe TO SELF: using a json file connect the grid specifications with both cpp and py files 
nx, ny = 64, 64
Lx, Ly = 1.0, 1.0
x = np.linspace(-Lx/2, Lx/2, nx)
y = np.linspace(-Ly/2, Ly/2, ny)
X, Y = np.meshgrid(x, y)

# heatmap of energy E
plt.figure(figsize=(10, 8))
heatmap = plt.imshow(energy_data, 
                    extent=[-Lx/2, Lx/2, -Ly/2, Ly/2],
                    cmap='viridis',
                    origin='lower')


cbar = plt.colorbar(heatmap, label='Energy (E)')
cbar.set_label('Energy (E)', rotation=270, labelpad=20)

# Axis labels and title
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('2D Gaussian Energy Distribution')

# Save and show
# plt.savefig('2d_gaussian.png', dpi=300, bbox_inches='tight')
plt.show()




# ###################### 3D plot work in progress ######################
# with h5py.File('energy_out.h5', 'r') as f:
#     dataset = f['Energy'][:]
#     energy_data = dataset[:, :, 0]  # First component is energy

# print("Data verification:")
# print(f"Min: {np.min(energy_data):.4f}, Max: {np.max(energy_data):.4f}")
# print(f"Center value: {energy_data[32,32]:.4f} (should be ~1.0)")

# nx, ny = 64, 64
# x = np.linspace(-0.5, 0.5, nx)
# y = np.linspace(-0.5, 0.5, ny)
# X, Y = np.meshgrid(x, y)

# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_surface(X, Y, energy_data, cmap='viridis', 
#                       rstride=1, cstride=1, antialiased=True)

# ax.set_xlabel('X Position')
# ax.set_ylabel('Y Position')
# ax.set_zlabel('Energy Density')
# plt.title('3D Gaussian Energy Distribution')
# fig.colorbar(surf, shrink=0.5, aspect=10)

# #plt.savefig('gaussian_3d.png', dpi=300, bbox_inches='tight')
# plt.show()
