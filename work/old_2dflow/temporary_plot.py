# energy_analysis.py
import h5py
import numpy as np
import matplotlib.pyplot as plt
import json

def analyze_simulation_results(hdf5_path, json_path):
    """Analyzes and visualizes simulation results with error handling."""
    try:
        # Load simulation parameters with validation
        with open(json_path) as f:
            config = json.load(f)
            
        required_keys = {
            'time_settings': ['t_start', 't_end', 'dt'],
            'grid': ['nx', 'ny'],
            'physical_size': ['Lx', 'Ly'],
            'initial_conditions': ['amplitude', 'sigma']
        }
        
        # Validate JSON structure
        for section, keys in required_keys.items():
            if section not in config:
                raise ValueError(f"Missing section: {section}")
            for key in keys:
                if key not in config[section]:
                    raise ValueError(f"Missing key: {section}.{key}")

        params = {
            'E0': config['initial_conditions']['amplitude'],
            'sigma': config['initial_conditions']['sigma'],
            't_start': config['time_settings']['t_start'],
            'nx': config['grid']['nx'],
            'ny': config['grid']['ny']
        }

        # Load and process HDF5 data
        with h5py.File(hdf5_path, 'r') as f:
            # Extract and sort time points
            datasets = sorted(f.keys(), key=lambda x: float(x.split('_t_')[1]))
            times = [float(name.split('_t_')[1]) for name in datasets]
            
            # Initialize storage
            avg_energy = np.zeros(len(datasets))
            center_energy = np.zeros(len(datasets))
            
            # Process each timestep
            for i, name in enumerate(datasets):
                data = np.array(f[name])
                
                # Handle 3D data (components, y, x)
                if data.ndim == 3:
                    e_data = data[0]  # Extract energy component
                else:
                    e_data = data.squeeze()
                
                # Calculate metrics
                avg_energy[i] = np.mean(e_data)
                
                # Get center value (handles 1D/2D/3D)
                if e_data.ndim == 2:
                    cy, cx = e_data.shape[0]//2, e_data.shape[1]//2
                    center_energy[i] = e_data[cy, cx]
                else:
                    center_energy[i] = e_data[len(e_data)//2]

        # Calculate analytical solution
        analytical = params['E0'] * (np.array(times)/params['t_start'])**(-4/3)

        # Create visualization
        plt.figure(figsize=(12, 7))
        
        # Plot numerical results
        plt.plot(times, avg_energy, 'b-', lw=2, label='Numerical (Average)')
        plt.plot(times, center_energy, 'r--', lw=2, label='Numerical (Center)')
        
        # Plot analytical solution
        plt.plot(times, analytical, 'k:', lw=3, label=r'Analytical $E_0(t/t_0)^{-4/3}$')
        
        # Format plot
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Energy', fontsize=14)
        plt.title('Energy Evolution Analysis', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=12)
        
        # Add parameter box
        text_content = [
            f'Initial Amplitude: {params["E0"]}',
            f'Start Time (t0): {params["t_start"]}',
            f'sigma: {params["sigma"]}',
            f'Grid: {params["nx"]}x{params["ny"]}'
        ]
        plt.gca().text(
            0.95, 0.95, 
            '\n'.join(text_content),
            transform=plt.gca().transAxes,
            fontsize=11,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(facecolor='white', alpha=0.8)
        )
        
        plt.tight_layout()
        plt.show()

    except FileNotFoundError as e:
        print(f"File error: {str(e)}")
    except json.JSONDecodeError as e:
        print(f"JSON decode error: {str(e)}")
    except KeyError as e:
        print(f"Missing key in config: {str(e)}")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    analyze_simulation_results('energy_out.h5', '2dflow_input.json')
