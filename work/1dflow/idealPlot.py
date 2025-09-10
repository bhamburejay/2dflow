import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_hdf5(filename, field="E", timestep=None):
    with h5py.File(filename, "r") as f:
        datasets = list(f.keys())
        field_map = {"E": 0, "M": 1, "e": 2, "ux": 3, "p": 4, "beta": 5, "cs2": 6}
        if isinstance(field, str):
            idx = field_map.get(field, 0)
        else:
            idx = int(field)
        initial_names = ["initialdata", "initialdatain"]
        for initial_name in initial_names:
            if initial_name in datasets:
                data = f[initial_name]
                if data.ndim == 3:
                    arr = data[0, :, :]
                elif data.ndim == 2:
                    arr = data[:, :]
                else:
                    arr = None
                if arr is not None and idx < arr.shape[1]:
                    y_init = arr[:, idx]
                    plt.plot(y_init, label=f"initial {field}", linestyle="--")
                break
        if "finaldata" in datasets:
            data = f["finaldata"]
            if data.ndim == 3:
                arr = data[-1, :, :]
            elif data.ndim == 2:
                arr = data[:, :]
            else:
                arr = None
            if arr is not None and idx < arr.shape[1]:
                y_final = arr[:, idx]
                plt.plot(y_final, label=f"final {field}")
        if ("initialdata" not in datasets) and ("finaldata" not in datasets):
            for name in ["solution"]:
                if name in datasets:
                    data = f[name]
                    if data.ndim == 3:
                        if timestep is None:
                            timestep = -1
                        arr = data[timestep, :, :]
                    elif data.ndim == 2:
                        arr = data[:, :]
                    else:
                        arr = None
                    if arr is not None and idx < arr.shape[1]:
                        y = arr[:, idx]
                        plt.plot(y, label=f"{field}")
        plt.xlabel("Grid index")
        plt.ylabel(field)
        plt.title(f"{field} (initial vs final)")
        plt.legend()
        plt.tight_layout()
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot 1D ideal hydro output from HDF5 file.")
    parser.add_argument("filename", help="HDF5 file to plot")
    parser.add_argument("--field", default="E", help="Field to plot (E, M, e, ux, p, beta, cs2)")
    parser.add_argument("--timestep", type=int, default=None, help="Timestep to plot (for time-dependent data)")
    args = parser.parse_args()
    plot_hdf5(args.filename, args.field, args.timestep)

if __name__ == "__main__":
    main()