import h5py
import numpy as np
import matplotlib.pyplot as plt

def plot_energy_vs_time(hdf5_path):
    with h5py.File(hdf5_path, 'r') as f:
        # Get sorted datasets by time
        datasets = sorted(f.keys(), key=lambda x: float(x.split('_t_')[1]))
        
        # Initialize arrays
        times = np.zeros(len(datasets))
        avg_energy = np.zeros(len(datasets))
        center_energy = np.zeros(len(datasets))
        
        for i, name in enumerate(datasets):
            times[i] = float(name.split('_t_')[1])
            data = f[name][:].squeeze()
            
            # Handle 3D data (time, y, x)
            if data.ndim == 3:
                # Take first time slice and find 2D center
                slice_data = data[0]
                ny, nx = slice_data.shape
                avg_energy[i] = np.mean(slice_data)
                center_energy[i] = slice_data[ny//2, nx//2]
            elif data.ndim == 2:
                # Standard 2D case
                ny, nx = data.shape
                avg_energy[i] = np.mean(data)
                center_energy[i] = data[ny//2, nx//2]
            else:
                # Handle 1D case
                avg_energy[i] = np.mean(data)
                center_energy[i] = data[len(data)//2]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(times, avg_energy, 'b-', label='Average Energy')
    plt.plot(times, center_energy, 'r--', label='Center Energy')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.title('Energy vs Time')
    plt.legend()
    plt.grid(True)
    plt.show()

plot_energy_vs_time('energy_out.h5')
