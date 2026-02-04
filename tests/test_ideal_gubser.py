# Script to run and analyze Gubser flow verification test
#
# This script sets up initial conditions for Gubser flow in a 2D viscous hydrodynamics simulation,
# runs the simulation, and then analyzes and plots the results comparing numerical and analytic solutions.
#
# To run the code
#
# python tests/test_ideal_gubser.py run --run_command "./VischydroMain.exe"
#
# To plot results
#
# python tests/test_ideal_gubser.py plot
import copy
import numpy as np
import matplotlib.pyplot as plt
import h5py
import subprocess
import os
import sys

# Import the helper module
import vh as vh


def gubser_solution(tau, r, q=1.0, e0=1.0):
    """
    Analytic solution for Gubser flow.
    """
    # Calculate rho (flow rapidity)
    numerator = 2 * q**2 * tau * r
    denominator = 1 + q**2 * (tau**2 + r**2)

    # Avoid domain errors in arctanh
    arg = numerator / denominator
    arg = np.clip(arg, -1.0 + 1e-15, 1.0 - 1e-15)
    rho = np.arctanh(arg)

    # Calculate energy density
    factor1 = e0 * (2 * q)**(8.0/3.0) / tau**(4.0/3.0)
    term2 = 1 + 2 * q**2 * (tau**2 + r**2) + q**4 * (tau**2 - r**2)**2
    e = factor1 / term2**(4.0/3.0)

    # Velocity components (u^tau, u^r)
    # u_tau = cosh(rho), u_r = sinh(rho)
    ur = np.sinh(rho)
    utau = np.cosh(rho)

    return e, utau, ur


<<<<<<< HEAD
def setup(q, e0, tau0):
    # Grid setup
    nx = 80
    ny = 80
=======
def run(run_name, q, e0, tau0, exe_path="./VischydroMain.exe"):

    # Grid setup - COMPROMISE: medium resolution
    nx = 60  # Between 40 and 80
    ny = 60  # Between 40 and 80
>>>>>>> 9a441c8 (gubser ideal test updated)
    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0
    vh.input_data["Vischydro/is_bjorken"] = True
    vh.input_data["Vischydro/use_ideal_step_only"] = True

    vh.input_data["VischydroMain/eta_by_s"] = 1.0 / (4.*np.pi)
    vh.input_data["VischydroMain/zeta_by_s"] = 0.5 / (4.*np.pi)
    vh.input_data["VischydroMain/initial_field_type"] = "charges"

    vh.input_data["VischydroMain/nx"] = nx
    vh.input_data["VischydroMain/ny"] = ny
    vh.input_data["VischydroMain/xmin"] = xmin
    vh.input_data["VischydroMain/xmax"] = xmax
    vh.input_data["VischydroMain/ymin"] = ymin
    vh.input_data["VischydroMain/ymax"] = ymax

    # Set time stepping - COMPROMISE: reasonable speed vs accuracy
    dt = 0.075  # Between 0.05 and 0.1
    t_end = 4.0  # Between 3.0 and 5.0
    vh.input_data["VischydroMain/t_start"] = tau0
    vh.input_data["VischydroMain/t_end"] = t_end
    vh.input_data["VischydroMain/dt_max"] = dt
    vh.input_data["VischydroMain/print_frequency"] = 15  # Between 10 and 20
    vh.input_data["VischydroMain/run_name"] = "gubser_test"


def run(q, e0, tau0, exe_path="./VischydroMain.exe"):
    setup(q, e0, tau0)

    # Generate initial condition
    X, Y, dx, dy = vh.xygrid(vh.input_data)

    R = np.sqrt(X**2 + Y**2)
    phi = np.arctan2(Y, X)

    print("Generating initial conditions...")

    e_init, utau_init, ur_init = gubser_solution(tau0, R, q, e0)
    ux_init = ur_init * np.cos(phi)
    uy_init = ur_init * np.sin(phi)

    # # The C++ code will:
    # # 1. Load the file.
    # # 2. Re-calculate primitive variables FROM the Conserved variables in the file.

    p_init = e_init / 3.0

    # Conserved variables
    E_cons = (e_init + p_init) * utau_init**2 - p_init
    Mx_cons = (e_init + p_init) * utau_init * ux_init
    My_cons = (e_init + p_init) * utau_init * uy_init

    initial_grid = vh.initialdata_grid(vh.input_data)
    initial_grid[:, :, 0] = E_cons
    initial_grid[:, :, 1] = Mx_cons
    initial_grid[:, :, 2] = My_cons
    initial_grid[:, :, 3] = e_init  # initial guess for root finder
    # these are not needed since we are using conserved initial conditions
    # initial_grid[:,:,4] = ux_init
    # initial_grid[:,:,5] = uy_init

    # print(f"Running simulation with {exe_path}...")
    vh.runcode(initial_grid, vh.input_data, runcommand=exe_path)


def analyze_results(q, e0, tau0):
    setup(q, e0, tau0)
    run_name = vh.input_data["VischydroMain/run_name"]

    # Load grid time file
    try:
        times = np.loadtxt(f"{run_name}_grid_t.txt")
        # Ensure times is iterable
        if times.ndim == 0:
            times = np.array([times])
    except:
        print("Could not load time file.")
        return

    # Load HDF5 grid
    h5_file = f"{run_name}_grid.h5"
    if not os.path.exists(h5_file):
        print("No output file found.")
        return

    with h5py.File(h5_file, 'r') as f:
        # Get coordinates
        # Assuming uniform grid, we can reconstruct or read
        # keys usually are 'coordinates', 'solution_step_N'

        # Plotting
        plt.figure(figsize=(10, 6))

        # Select a few steps to plot (e.g., all steps available)
        # HDF5 dataset "solution": shape (N_steps, ny, nx, ndof)

        N_steps = f['solution'].shape[0]
        # step_indices = np.linspace(0, N_steps-1, 5, dtype=int)
        step_indices = range(N_steps)

        print(times)
        print(times.shape)
        for idx in step_indices:
            time_val = times[idx, 0]

            # group_name = f"solution_step_{step_num}"
            # Data is in dataset 'solution' at index idx
            data = f['solution'][idx]  # Shape (ny, nx, ndof) usually

            # Extract center slice (y=0 approx)
            ny, nx, ndof = data.shape
            mid_y = ny // 2

            xs = f['coordinates'][mid_y, :, 0]
            ys = f['coordinates'][mid_y, :, 1]

            r_vals = np.sqrt(xs**2 + ys**2)

            # Calculate T00 = E
            # Analytic comparison

            e_ana, ut_ana, ur_ana = gubser_solution(time_val, r_vals, q, e0)

            # Plot tau * T00 approx = tau * E
            # Actually user asked for tau * T00
            # E in code is T00.

            E_num = data[mid_y, :, 0]

            label_txt = f"t={time_val:.2f}"
            print(f"Plotting for time {time_val:.2f}")
            plt.plot(xs, time_val * E_num, '.', label=f"Sim {label_txt}")

            # Analytic E
            # E = (e+p) u0^2 - p = (4/3 e) u0^2 - e/3 = e (4 u0^2 - 1)/3
            E_ana = e_ana * (4 * ut_ana**2 - 1.0) / 3.0
            plt.plot(xs, time_val * E_ana, '-',
                     alpha=0.5, label=f"Ana {label_txt}")

        plt.xlabel('r')
        plt.ylabel(r'$\tau T^{00}$')
        plt.title(f'Gubser Flow - Conserved Init Test')
        plt.legend()
        plt.grid(True)
        plt.savefig(f"gubser_ideal_test.png")
        print(f"Plot saved to gubser_ideal_test.png")


if __name__ == "__main__":
    import argparse
    # Two modes : gubser_test0 run --run_command <command> or gubser_test0 plot
    parser = argparse.ArgumentParser(
        description="Gubser Flow Verification Test", usage="%(prog)s {run|plot} [options]")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Sub-command: run --run_command <command>
    run_parser = subparsers.add_parser(
        "run", help="Run the Gubser flow verification test")

    run_parser.add_argument("--run_command", default="./VischydroMain.exe",
                            help="Command to run the VischydroMain executable")

    # Sub-command: plot
    plot_parser = subparsers.add_parser(
        "plot", help="Plot results from Gubser flow verification test")

    args = parser.parse_args()

    # Setup parameters
    q = 1.0
    e0 = 1.0
    tau0 = 0.6

    if args.command == "run":
        run(q, e0, tau0, exe_path=args.run_command)
    elif args.command == "plot":
        analyze_results(q, e0, tau0)
