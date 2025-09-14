import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

# Usage: python plot_gauss2d.py <output_file.h5> [field]
# Default field is 'E'

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Plot 2D hydrodynamics output from HDF5 file.')
    parser.add_argument('filename', help='HDF5 output file (e.g. *_initial.h5 or *_final.h5)')
    parser.add_argument('--field', default='E', help='Field to plot (default: E)')
    parser.add_argument('--step', type=int, default=None, help='Timestep to plot (if time series)')
    args = parser.parse_args()

    with h5py.File(args.filename, 'r') as f:
        # Try to find the solution dataset
        if 'solution' in f:
            data = f['solution'][:]
        elif args.field in f:
            data = f[args.field][:]
        else:
            # Try to find a group with the field
            for key in f.keys():
                if args.field in f[key]:
                    data = f[key][args.field][:]
                    break
            else:
                print(f"Field '{args.field}' not found in file.")
                return

        # If time series, select step
        if data.ndim == 3:
            if args.step is None:
                print(f"Data has {data.shape[0]} timesteps. Use --step to select.")
                step = data.shape[0] - 1
                print(f"Plotting last step: {step}")
            else:
                step = args.step
            arr = data[step]
        elif data.ndim == 2:
            arr = data
        else:
            print(f"Unexpected data shape: {data.shape}")
            return

        # Plot
        plt.figure(figsize=(8,6))
        plt.imshow(arr, origin='lower', aspect='auto', cmap='viridis')
        plt.colorbar(label=args.field)
        plt.title(f"{args.field} from {args.filename}")
        plt.xlabel('X index')
        plt.ylabel('Y index')
        plt.tight_layout()
        plt.show()

if __name__ == '__main__':
    main()
