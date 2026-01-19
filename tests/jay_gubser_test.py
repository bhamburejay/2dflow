import numpy as np
import matplotlib.pyplot as plt
import h5py
import subprocess
import os
import sys

# Import the helper module
import vischydro as vh

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

import copy

def run_verification():
    # Setup parameters
    q = 1.0
    e0 = 1.0
    tau0 = 1.0
    
    # Grid setup
    nx = 200
    ny = 3  # Small y dimension for 1D-like simulation (assuming check x-axis)
    xmin, xmax = -4.0, 4.0
    ymin, ymax = -0.1, 0.1
    
    # Configure input parameters
    # Must deepcopy and wrap to support runcode's slash access AND allow new keys
    input_data_raw = copy.deepcopy(vh._input_data)
    
    # Add the new key to the raw dictionary structure
    # Now passing directly in Vischydro config section
    input_data_raw["Vischydro"]["init_from_conserved_charged"] = True
    input_data_raw["VischydroMain"]["eta_by_s"] = 0.0 # Ideal hydro for Gubser match
    input_data_raw["VischydroMain"]["zeta_by_s"] = 0.0
    
    # Wrap in FixedData
    input_data = vh.FixedData(input_data_raw)
    
    # Set values via FixedData interface (or raw)
    input_data["Vischydro/nx"] = nx
    input_data["Vischydro/ny"] = ny
    input_data["Vischydro/xmin"] = xmin
    input_data["Vischydro/xmax"] = xmax
    input_data["Vischydro/ymin"] = ymin
    input_data["Vischydro/ymax"] = ymax
    input_data["Vischydro/is_bjorken"] = True
    
    # Set time stepping
    dt = 0.05
    t_end = 3.0
    input_data["VischydroMain/t_start"] = tau0
    input_data["VischydroMain/t_end"] = t_end
    input_data["VischydroMain/dt_max"] = dt
    input_data["VischydroMain/print_frequency"] = 10 
    input_data["VischydroMain/run_name"] = "gubser_test"


    # Generate initial condition
    X, Y = np.meshgrid(
        np.linspace(xmin, xmax, nx),
        np.linspace(ymin, ymax, ny)
    )
    R = np.sqrt(X**2 + Y**2)
    # For 1D slice along X, we can treat X as r if Y~0. 
    # But Gubser is cylindrically symmetric.
    # Let's simple initialize using r = |x|
    
    e_init, utau_init, ur_init = gubser_solution(tau0, np.abs(X), q, e0)
    
    # Map ur to ux, uy
    # u^x = u^r * (x/r)
    # Handle r=0
    with np.errstate(invalid='ignore'):
        cos_phi = X / R
        sin_phi = Y / R
    cos_phi[R==0] = 1.0
    sin_phi[R==0] = 0.0
    
    ux_init = ur_init * cos_phi
    uy_init = ur_init * sin_phi
    
    # ndof=9 assumed: [E, M1, M2, e, u1, u2, p, beta, cs2] ?? 
    # vhnode uses: e, u[0], u[1]...
    # The python helper vh.runcode usually expects a grid of primitives or full nodes?
    # Looking at vischydro.py, it likely initializes something.
    # Let's create the 'initialdata' array.
    # We need to fill it ourselves or let the C++ code fill it?
    # The C++ code loads from HDF5.
    
    # We will construct the primitive variable array [e, ux, uy, 0...]
    # and let the Code calculate Conserved variables, 
    # BUT since we set init_from_conserved_charged = True,
    # The C++ code will:
    # 1. Load the file.
    # 2. Re-calculate Primitive variables FROM the Conserved variables in the file.
    # Wait. If the file contains Primitive variables, and we built Conserved from them...
    # VischydroNode.cpp: fill() calculates E, M from e, u.
    # If we save that to file, the file has everything consistent.
    # Then init_from_conserved_charged means:
    # "Read E, M from file. Discard e, u. Re-solve e, u from E, M."
    # This validates the root finder.
    
    # We need to constructing the full node data structure to save to H5.
    # We can use the VischydroNode structure approximation.
    # Or just use the 'vh' helper if it has one.
    
    ndof = 9
    initial_grid = np.zeros((ny, nx, ndof))
    
    # Assume layout: 
    # 0: E, 1: Mx, 2: My, 3: e, 4: ux, 5: uy, 6: p, 7: beta, 8: cs2
    # This might vary based on struct layout, but usually:
    # struct VischydroNode { E, M[2], e, u[2], p, beta, cs2 }
    
    # Simple EOS P = e/3
    p_init = e_init / 3.0
    
    # Conserved variables
    # E = (e+p) u0^2 - p
    # M = (e+p) u0 u
    
    E_cons = (e_init + p_init) * utau_init**2 - p_init
    Mx_cons = (e_init + p_init) * utau_init * ux_init
    My_cons = (e_init + p_init) * utau_init * uy_init
    
    initial_grid[:,:,0] = E_cons
    initial_grid[:,:,1] = Mx_cons
    initial_grid[:,:,2] = My_cons
    initial_grid[:,:,3] = e_init
    initial_grid[:,:,4] = ux_init
    initial_grid[:,:,5] = uy_init
    initial_grid[:,:,6] = p_init
    # Fill others with dummy or approximation
    initial_grid[:,:,7] = 0.0 # beta
    initial_grid[:,:,8] = 1.0/3.0 # cs2

    # Run Simulation
    # Check if executable exists
    exe_path = "../build/tests/VischydroMain.exe"
    if not os.path.exists(exe_path):
        print(f"Executable not found at {exe_path}")
        # Try finding in current or standard build loc
        exe_path = "./VischydroMain.exe"
    
    print(f"Running simulation with {exe_path}...")
    vh.runcode(initial_grid, input_data, runcommand=exe_path)
    
    # Post-process results
    analyze_results(input_data["VischydroMain"]["run_name"], q, e0)

def analyze_results(run_name, q, e0):
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

        for idx in step_indices:
            # Get valid time from time file if possible, or infer
            if idx < len(times):
                step_num = int(times[idx, 0])
                time_val = times[idx, 1]
            else:
                print(f"Index {idx} out of bounds for times array.")
                continue

            # group_name = f"solution_step_{step_num}"
            # Data is in dataset 'solution' at index idx
            data = f['solution'][idx] # Shape (ny, nx, ndof) usually
            
            # Extract center slice (y=0 approx)
            ny, nx, ndof = data.shape
            mid_y = ny // 2
            
            # Get x coordinates (assuming standard grid)
            # We can infer from data size and hardcoded bounds for now
            # or read 'coordinates' from file if available
            xs = np.linspace(-4.0, 4.0, nx)
            
            e_num = data[mid_y, :, 3] # Primitive e
            utau_num = data[mid_y, :, 0] # Wait, solution usually stores VischydroNode layout
            # which is E, M, e, u ...
            # 0: E, 1: M0, 2: M1, 3:e, 4:u0, 5:u1
            
            # Calculate T00 = E
            # Analytic comparison
            
            e_ana, ut_ana, ur_ana = gubser_solution(time_val, np.abs(xs), q, e0)
            
            # Plot tau * T00 approx = tau * E
            # Actually user asked for tau * T00
            # E in code is T00.
            
            E_num = data[mid_y, :, 0]
            
            label_txt = f"t={time_val:.2f}"
            plt.plot(xs, time_val * E_num, '.', label=f"Sim {label_txt}")
            
            # Analytic E
            # E = (e+p) u0^2 - p = (4/3 e) u0^2 - e/3 = e (4 u0^2 - 1)/3
            E_ana = e_ana * (4 * ut_ana**2 - 1.0) / 3.0
            plt.plot(xs, time_val * E_ana, '-', alpha=0.5, label=f"Ana {label_txt}")

        plt.xlabel('r')
        plt.ylabel(r'$\tau T^{00}$')
        plt.title(f'Gubser Flow - Conserved Init Test')
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{run_name}_comparison.png")
        print(f"Plot saved to {run_name}_comparison.png")

if __name__ == "__main__":
    run_verification()
