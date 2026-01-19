import h5py
import matplotlib.pyplot as plt
import numpy as np

def plot_viscous_tensor():
    with h5py.File('simple_bj2d_example_grid.h5', 'r') as f:
        # Shape: (t, y, x, components)
        # Components: 4 -> xx, 5 -> xy, 7 -> yy
        pi = f['viscous_tensor'][:]
        coords = f['coordinates'][:]
        
        # Last time step
        pi_last = pi[-1, :, :, :]
        
        # Plot Pi_xx
        plt.figure(figsize=(10, 8))
        plt.imshow(pi_last[:, :, 4], origin='lower', cmap='RdBu_r')
        plt.colorbar(label=r'$\pi^{xx}$')
        plt.title(r'$\pi^{xx}$ at final time')
        plt.savefig('viscous_tensor_pixx.png')
        print("Saved viscous_tensor_pixx.png")

        # Plot Pi_xy
        plt.figure(figsize=(10, 8))
        plt.imshow(pi_last[:, :, 5], origin='lower', cmap='RdBu_r')
        plt.colorbar(label=r'$\pi^{xy}$')
        plt.title(r'$\pi^{xy}$ at final time')
        plt.savefig('viscous_tensor_pixy.png')
        print("Saved viscous_tensor_pixy.png")

if __name__ == "__main__":
    plot_viscous_tensor()
