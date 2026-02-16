# Script to run and analyze Bjorken flow verification test
#
# This script sets up initial conditions for Bjorken flow in a 1D boost-invariant
# hydrodynamics simulation, runs the simulation, and then analyzes and plots the results
# comparing numerical and analytic solutions.
#
# To run the code
#
# python tests/test_bjorken_flow.py run --run_command "./VischydroMain.exe"
#
# To plot results
#
# python tests/test_bjorken_flow.py plot
import copy
import numpy as np
import matplotlib.pyplot as plt
import h5py
import subprocess
import os
import sys

# Import the helper module
import vischydro as vh


def bjorken_solution(tau, e0=1.0, tau0=0.6):
    """
    Analytic solution for Bjorken flow (1D boost-invariant expansion).
    Energy density scales as τ^{-4/3}.
    """
    return e0 * (tau0 / tau)**(4.0/3.0)


def run(run_name, e0, tau0, exe_path="./VischydroMain.exe"):

    # Grid setup - 1D in practice (radial symmetry)
    nx = 60
    ny = 60
    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0

    vh.input_data["VischydroMain/eta_by_s"] = 0.0
    vh.input_data["VischydroMain/zeta_by_s"] = 0.0
    vh.input_data["VischydroMain/initial_field_type"] = "charges"  # Use conserved variables like Gubser

    vh.input_data["VischydroMain/nx"] = nx
    vh.input_data["VischydroMain/ny"] = ny
    vh.input_data["VischydroMain/xmin"] = xmin
    vh.input_data["VischydroMain/xmax"] = xmax
    vh.input_data["VischydroMain/ymin"] = ymin
    vh.input_data["VischydroMain/ymax"] = ymax
    vh.input_data["Vischydro/is_bjorken"] = True

    # Set time stepping
    dt = 0.075
    t_end = 4.0
    vh.input_data["VischydroMain/t_start"] = tau0
    vh.input_data["VischydroMain/t_end"] = t_end
    vh.input_data["VischydroMain/dt_max"] = dt
    vh.input_data["VischydroMain/print_frequency"] = 15
    vh.input_data["VischydroMain/run_name"] = run_name

    # Generate initial condition - uniform energy density for Bjorken
    X, Y, dx, dy = vh.xygrid(vh.input_data)

    print("Generating initial conditions for Bjorken flow...")

    # For Bjorken flow, we need conserved variables
    # E = (ε + p)γ² - p, but for initial rest frame, γ=1, so E = ε + p
    e_init = e0 * np.ones_like(X)  # energy density
    p_init = e_init / 3.0  # pressure for ideal gas

    # Conserved variables for Bjorken flow (boost-invariant)
    # In lab frame initially: E = (ε + p) - p = ε
    E_cons = e_init  # T^{00} = ε initially
    Mx_cons = np.zeros_like(X)  # No transverse flow
    My_cons = np.zeros_like(X)

    initial_grid = vh.initialdata_grid(vh.input_data)
    initial_grid[:, :, 0] = E_cons
    initial_grid[:, :, 1] = Mx_cons
    initial_grid[:, :, 2] = My_cons
    initial_grid[:, :, 3] = e_init  # initial guess for root finder

    # print(f"Running simulation with {exe_path}...")
    vh.runcode(initial_grid, vh.input_data, runcommand=exe_path)


def analyze_results(run_name, e0, tau0):
    # Load grid time file
    try:
        times = np.loadtxt(f"{run_name}_grid_t.txt")
        # Ensure times is iterable
        if times.ndim == 0:
            times = np.array([times])
        elif times.ndim == 1 and len(times) > 1:
            # If it's a 1D array with multiple values, treat as single time step
            times = np.array([times[0]])  # Take first value as time
    except:
        print("Could not load time file.")
        return

    # Load HDF5 grid
    h5_file = f"{run_name}_grid.h5"
    if not os.path.exists(h5_file):
        print("No output file found.")
        return

    with h5py.File(h5_file, 'r') as f:
        # Standard Bjorken plot: ε vs τ on log-log scale
        plt.figure(figsize=(10, 6))

        N_steps = f['solution'].shape[0]

        print("Time steps:", times)

        # Collect data for plotting
        time_vals = []
        energy_vals = []
        analytic_vals = []

        for idx in range(N_steps):
            time_val = times[idx] if times.ndim == 1 else times[idx, 0]

            data = f['solution'][idx]

            # Extract center slice and average over transverse directions
            ny, nx, ndof = data.shape
            mid_y = ny // 2

            # For Bjorken flow, average energy density over transverse plane
            e_num = np.mean(data[mid_y, :, 0])  # T⁰⁰ ≈ ε for Bjorken

            # Analytic solution
            e_ana = bjorken_solution(time_val, e0, tau0)

            time_vals.append(time_val)
            energy_vals.append(e_num)
            analytic_vals.append(e_ana)

            print(f"Time {time_val:.2f}: ε_num={e_num:.4f}, ε_ana={e_ana:.4f}")

        # Plot on log-log scale (standard for Bjorken)
        plt.loglog(time_vals, energy_vals, 'ro-', markersize=8, linewidth=2,
                   label='Numerical Simulation')
        plt.loglog(time_vals, analytic_vals, 'b--', linewidth=3,
                   label=r'Analytic: $\epsilon \propto \tau^{-4/3}$')

        plt.xlabel(r'Proper Time $\tau$')
        plt.ylabel(r'Energy Density $\epsilon$')
        plt.title('Bjorken Flow: Energy Density Evolution')
        plt.legend()
        plt.grid(True, which='both', alpha=0.3)

        # Add slope reference line (τ^{-4/3})
        tau_ref = np.array([0.5, 5.0])
        slope_ref = e0 * (tau0 / tau_ref)**(4.0/3.0)
        plt.loglog(tau_ref, slope_ref, 'k:', alpha=0.5,
                   label=r'$\tau^{-4/3}$ reference')

        plt.tight_layout()
        plt.savefig(f"{run_name}_comparison.png", dpi=150, bbox_inches='tight')
        print(f"Plot saved to {run_name}_comparison.png")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Bjorken Flow Verification Test", usage="%(prog)s {run|plot} [options]")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Sub-command: run
    run_parser = subparsers.add_parser(
        "run", help="Run the Bjorken flow verification test")

    run_parser.add_argument("--run_command", required=True, type=str,
                            default="./VischydroMain.exe", help="Command to run the VischydroMain executable")

    # Sub-command: plot
    plot_parser = subparsers.add_parser(
        "plot", help="Plot results from Bjorken flow verification test")

    args = parser.parse_args()

    # Setup parameters
    run_name = "bjorken_test"
    e0 = 1.0
    tau0 = 0.6

    if args.command == "run":
        run(run_name, e0, tau0, exe_path=args.run_command)
    elif args.command == "plot":
        analyze_results(run_name, e0, tau0)