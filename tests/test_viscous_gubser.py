# Script to run and analyze viscous Gubser flow verification test
#
# This is intentionally structured like tests/test_ideal_gubser.py:
# - `run` builds initial conditions and launches VischydroMain
# - `plot` compares simulation against analytic curves
#
# Paper-oriented defaults (Gubser 2010):
#   1/q = 4.3 fm, eta/s = 0.134, hat{T}_0 = 5.55, tau0 = 0.5 fm/c
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import shutil
from scipy.special import hyp2f1

import vischydro as vh


OUTDIR = os.path.abspath(os.path.dirname(__file__))

DEFAULT_Q = 1.0 / 4.3
DEFAULT_ETA_BY_S = 0.134
DEFAULT_HAT_T0 = 5.55
DEFAULT_TAU0 = 0.5
DEFAULT_T_END = 9.0
DEFAULT_F_STAR = 11.0


def print_parameter_summary(title, params):
    print(f"\n{title}")
    for key, value in params.items():
        print(f"  {key}: {value}")


def resolve_executable(exe_path):
    if os.path.isabs(exe_path) or os.path.exists(exe_path):
        return exe_path
    if shutil.which(exe_path) is not None:
        return exe_path

    repo_root = os.path.abspath(os.path.join(OUTDIR, ".."))
    fallback_candidates = [
        os.path.join(repo_root, "build", "tests", "VischydroMain.exe"),
        os.path.join(OUTDIR, "VischydroMain.exe"),
    ]
    for candidate in fallback_candidates:
        if os.path.exists(candidate):
            print(f"Using fallback executable: {candidate}")
            return candidate

    raise FileNotFoundError(
        f"Could not find executable '{exe_path}'. "
        "Tried provided path and fallbacks: build/tests/VischydroMain.exe, tests/VischydroMain.exe"
    )


def viscous_gubser_solution(
    tau,
    r,
    q=DEFAULT_Q,
    hatT0=DEFAULT_HAT_T0,
    eta_by_s=DEFAULT_ETA_BY_S,
    f_star=DEFAULT_F_STAR,
):
    """
    Viscous Gubser solution using Eq.(26) in Gubser (2010).
    Returns: e, utau, ur, hatT.
    """
    tau = float(tau)
    r = np.asarray(r, dtype=float)

    numerator = 2.0 * q * q * tau * r
    denominator = 1.0 + q * q * (tau * tau + r * r)
    v_perp = numerator / (denominator + 1e-30)
    v_perp = np.clip(v_perp, -1.0 + 1e-15, 1.0 - 1e-15)
    rho = np.arctanh(v_perp)
    utau = np.cosh(rho)
    ur = np.sinh(rho)

    g = (1.0 + q * q * r * r - q * q * tau * tau) / (2.0 * q * tau + 1e-30)
    z = -g * g
    hyp = hyp2f1(0.5, 1.0 / 6.0, 1.5, z)

    sigma0 = 4.0 / 3.0 * f_star**0.25
    H0 = eta_by_s * sigma0
    hatT = (
        hatT0 / (1.0 + g * g) ** (1.0 / 3.0)
        + (H0 * g / np.sqrt(1.0 + g * g))
        * (1.0 - (1.0 + g * g) ** (1.0 / 6.0) * hyp)
    )

    # e = hat{T}^4 / tau^4; invalid region has hatT <= 0 in first-order hydro.
    e = np.where(hatT > 0.0, (hatT**4) / (tau**4 + 1e-30), np.nan)
    return e, utau, ur, hatT


def ideal_gubser_solution(tau, r, q, e0):
    tau = float(tau)
    r = np.asarray(r, dtype=float)

    numerator = 2.0 * q * q * tau * r
    denominator = 1.0 + q * q * (tau * tau + r * r)
    arg = numerator / (denominator + 1e-30)
    arg = np.clip(arg, -1.0 + 1e-15, 1.0 - 1e-15)
    rho = np.arctanh(arg)

    factor1 = e0 * (2.0 * q) ** (8.0 / 3.0) / (tau ** (4.0 / 3.0) + 1e-30)
    term2 = 1.0 + 2.0 * q * q * (tau * tau + r * r) + q**4 * (tau * tau - r * r) ** 2
    e = factor1 / (term2 ** (4.0 / 3.0) + 1e-30)
    utau = np.cosh(rho)
    ur = np.sinh(rho)
    return e, utau, ur


def run(
    run_name,
    q,
    hatT0,
    tau0,
    eta_by_s=DEFAULT_ETA_BY_S,
    f_star=DEFAULT_F_STAR,
    t_end=DEFAULT_T_END,
    dt=0.075,
    print_frequency=1,
    nx=60,
    ny=60,
    xmin=-10.0,
    xmax=10.0,
    ymin=-10.0,
    ymax=10.0,
    exe_path="./VischydroMain.exe",
):
    run_params = {
        "q [1/fm]": q,
        "1/q [fm]": (1.0 / q) if q != 0.0 else np.inf,
        "eta/s": eta_by_s,
        "hatT0": hatT0,
        "f_star": f_star,
        "tau0 [fm/c]": tau0,
        "t_end [fm/c]": t_end,
        "dt_max": dt,
        "print_frequency": print_frequency,
        "nx": nx,
        "ny": ny,
        "xmin": xmin,
        "xmax": xmax,
        "ymin": ymin,
        "ymax": ymax,
        "run_name": run_name,
    }
    print_parameter_summary("Run Parameters", run_params)

    vh.input_data["VischydroMain/eta_by_s"] = eta_by_s
    vh.input_data["VischydroMain/zeta_by_s"] = 0.0
    vh.input_data["VischydroMain/initial_field_type"] = "charges"

    vh.input_data["VischydroMain/nx"] = nx
    vh.input_data["VischydroMain/ny"] = ny
    vh.input_data["VischydroMain/xmin"] = xmin
    vh.input_data["VischydroMain/xmax"] = xmax
    vh.input_data["VischydroMain/ymin"] = ymin
    vh.input_data["VischydroMain/ymax"] = ymax
    vh.input_data["Vischydro/is_bjorken"] = True

    vh.input_data["VischydroMain/t_start"] = tau0
    vh.input_data["VischydroMain/t_end"] = t_end
    vh.input_data["VischydroMain/dt_max"] = dt
    vh.input_data["VischydroMain/print_frequency"] = print_frequency
    vh.input_data["VischydroMain/run_name"] = os.path.join(OUTDIR, run_name)

    X, Y, dx, dy = vh.xygrid(vh.input_data)
    R = np.sqrt(X**2 + Y**2)
    phi = np.arctan2(Y, X)

    print("Generating viscous initial conditions...")
    e_init, utau_init, ur_init, hatT_init = viscous_gubser_solution(
        tau0, R, q=q, hatT0=hatT0, eta_by_s=eta_by_s, f_star=f_star
    )
    ux_init = ur_init * np.cos(phi)
    uy_init = ur_init * np.sin(phi)

    # Keep initialization robust where first-order hydro is not valid (hatT<=0).
    bad = ~np.isfinite(e_init)
    e_floor = 1e-12
    if np.any(bad):
        e_init = e_init.copy()
        utau_init = utau_init.copy()
        ux_init = ux_init.copy()
        uy_init = uy_init.copy()
        e_init[bad] = e_floor
        utau_init[bad] = 1.0
        ux_init[bad] = 0.0
        uy_init[bad] = 0.0
    e_init = np.maximum(e_init, e_floor)

    p_init = e_init / 3.0
    E_cons = (e_init + p_init) * utau_init**2 - p_init
    Mx_cons = (e_init + p_init) * utau_init * ux_init
    My_cons = (e_init + p_init) * utau_init * uy_init

    initial_grid = vh.initialdata_grid(vh.input_data)
    initial_grid[:, :, 0] = E_cons
    initial_grid[:, :, 1] = Mx_cons
    initial_grid[:, :, 2] = My_cons
    initial_grid[:, :, 3] = e_init  # root-finder initial guess

    exe_candidate = resolve_executable(exe_path)
    vh.runcode(initial_grid, vh.input_data, runcommand=exe_candidate)


def get_time_column(times):
    arr = np.asarray(times, dtype=float)
    if arr.ndim == 0:
        return np.array([float(arr)])
    if arr.ndim == 1:
        # one-row file with columns [time, step, save_index]
        if arr.size >= 3:
            return np.array([float(arr[0])])
        return arr
    return arr[:, 0]


def analyze_results(
    run_name,
    q,
    hatT0,
    eta_by_s=DEFAULT_ETA_BY_S,
    f_star=DEFAULT_F_STAR,
    max_curves=6,
    overlay_ideal=True,
    ideal_run_name="gubser_test",
):
    plot_params = {
        "q [1/fm]": q,
        "1/q [fm]": (1.0 / q) if q != 0.0 else np.inf,
        "eta/s": eta_by_s,
        "hatT0": hatT0,
        "f_star": f_star,
        "max_curves": max_curves,
        "overlay_ideal": overlay_ideal,
        "ideal_run_name": ideal_run_name,
        "run_name": run_name,
    }
    print_parameter_summary("Plot Parameters", plot_params)

    times_file = os.path.join(OUTDIR, f"{run_name}_grid_t.txt")
    try:
        times = np.loadtxt(times_file)
    except Exception:
        print(f"Could not load time file: {times_file}")
        return

    h5_file = os.path.join(OUTDIR, f"{run_name}_grid.h5")
    if not os.path.exists(h5_file):
        print(f"No output file found: {h5_file}")
        return

    time_vals = get_time_column(times)

    with h5py.File(h5_file, "r") as f:
        N_steps = f["solution"].shape[0]
        ny = f["solution"].shape[1]
        nx = f["solution"].shape[2]
        print(f"HDF5 solution shape: steps={N_steps}, ny={ny}, nx={nx}")

        if len(time_vals) != N_steps:
            nplot = min(len(time_vals), N_steps)
            print(
                f"Warning: times length ({len(time_vals)}) != N_steps ({N_steps}). "
                f"Using first {nplot} snapshots."
            )
        else:
            nplot = N_steps

        if nplot <= max_curves:
            step_indices = range(nplot)
        else:
            step_indices = np.unique(
                np.linspace(0, nplot - 1, num=max_curves, dtype=int)
            )
        mid_y = ny // 2
        xs = f["coordinates"][mid_y, :, 0]
        ys = f["coordinates"][mid_y, :, 1]
        r_vals = np.sqrt(xs**2 + ys**2)

        plt.figure(figsize=(10, 6))
        for idx in step_indices:
            tau = float(time_vals[idx])
            data = f["solution"][idx]
            E_num = data[mid_y, :, 0]

            e_ana, ut_ana, ur_ana, hatT_ana = viscous_gubser_solution(
                tau, r_vals, q=q, hatT0=hatT0, eta_by_s=eta_by_s, f_star=f_star
            )
            E_ana = e_ana * (4.0 * ut_ana**2 - 1.0) / 3.0
            valid = np.isfinite(E_ana)

            label_txt = rf"$\tau={tau:.2f}$"
            plt.plot(xs, tau * E_num, ".", label=f"Sim {label_txt}")
            plt.plot(
                xs[valid],
                tau * E_ana[valid],
                "-",
                alpha=0.6,
                label=f"Ana {label_txt}",
            )

        plt.xlabel("r")
        plt.ylabel(r"$\tau T^{00}$")
        plt.title("Viscous Gubser Flow - Conserved Init Test")
        plt.legend()
        plt.grid(True)
        out_png = os.path.join(OUTDIR, "gubser_viscous_test.png")
        plt.savefig(out_png)
        print(f"Plot saved to {out_png}")

    if not overlay_ideal:
        return

    # Optional overlay with ideal run if files are available.
    ideal_prefix = ideal_run_name or "gubser_test"
    ideal_h5_candidates = [
        os.path.join(OUTDIR, f"{ideal_prefix}_grid.h5"),
        os.path.join(os.path.abspath(os.path.join(OUTDIR, "..")), f"{ideal_prefix}_grid.h5"),
    ]
    ideal_t_candidates = [
        os.path.join(OUTDIR, f"{ideal_prefix}_grid_t.txt"),
        os.path.join(os.path.abspath(os.path.join(OUTDIR, "..")), f"{ideal_prefix}_grid_t.txt"),
    ]

    ideal_h5 = next((p for p in ideal_h5_candidates if os.path.exists(p)), None)
    ideal_t = next((p for p in ideal_t_candidates if os.path.exists(p)), None)
    if ideal_h5 is None or ideal_t is None:
        print("Ideal overlay skipped: ideal output files not found.")
        return

    try:
        ideal_times_raw = np.loadtxt(ideal_t)
        ideal_time_vals = get_time_column(ideal_times_raw)
    except Exception:
        print(f"Ideal overlay skipped: could not load {ideal_t}.")
        return

    with h5py.File(h5_file, "r") as fv, h5py.File(ideal_h5, "r") as fi:
        nv = fv["solution"].shape[0]
        ni = fi["solution"].shape[0]
        nv_use = min(nv, len(time_vals))
        ni_use = min(ni, len(ideal_time_vals))
        if nv_use == 0 or ni_use == 0:
            print("Ideal overlay skipped: one of the runs has no saved snapshots.")
            return

        tv = time_vals[:nv_use]
        ti = ideal_time_vals[:ni_use]

        # Pick a few viscous times and match nearest ideal times.
        nsel = min(max_curves, nv_use)
        v_indices = np.unique(np.linspace(0, nv_use - 1, num=nsel, dtype=int))

        ny_v = fv["solution"].shape[1]
        mid_y_v = ny_v // 2
        xs_v = fv["coordinates"][mid_y_v, :, 0]
        ys_v = fv["coordinates"][mid_y_v, :, 1]
        r_v = np.sqrt(xs_v**2 + ys_v**2)

        ny_i = fi["solution"].shape[1]
        mid_y_i = ny_i // 2
        xs_i = fi["coordinates"][mid_y_i, :, 0]
        ys_i = fi["coordinates"][mid_y_i, :, 1]
        r_i = np.sqrt(xs_i**2 + ys_i**2)

        # Estimate ideal e0 from first ideal snapshot at r~0.
        tau0_i = float(ti[0])
        i0 = 0
        e_center_i = float(fi["solution"][i0, mid_y_i, np.argmin(np.abs(r_i)), 3])
        denom0 = (2.0 * q) ** (8.0 / 3.0)
        term20 = 1.0 + 2.0 * q * q * (tau0_i * tau0_i) + q**4 * (tau0_i * tau0_i) ** 2
        e0_ideal = e_center_i * (tau0_i ** (4.0 / 3.0)) * (term20 ** (4.0 / 3.0)) / (denom0 + 1e-30)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
        ax_i, ax_v = axes
        for iv in v_indices:
            tau_v = float(tv[iv])
            ii = int(np.argmin(np.abs(ti - tau_v)))
            tau_i = float(ti[ii])

            Ev = fv["solution"][iv, mid_y_v, :, 0]
            Ei = fi["solution"][ii, mid_y_i, :, 0]

            # Ideal panel: simulation + ideal analytic line.
            e_ana_i, ut_ana_i, ur_ana_i = ideal_gubser_solution(tau_i, r_i, q=q, e0=e0_ideal)
            E_ana_i = e_ana_i * (4.0 * ut_ana_i**2 - 1.0) / 3.0
            ax_i.plot(r_i, tau_i * Ei, ".", label=rf"Sim $\tau={tau_i:.2f}$")
            ax_i.plot(r_i, tau_i * E_ana_i, "-", alpha=0.7, label=rf"Ana $\tau={tau_i:.2f}$")

            # Viscous panel: simulation + viscous analytic line.
            e_ana_v, ut_ana_v, ur_ana_v, hatT_ana_v = viscous_gubser_solution(
                tau_v, r_v, q=q, hatT0=hatT0, eta_by_s=eta_by_s, f_star=f_star
            )
            E_ana_v = e_ana_v * (4.0 * ut_ana_v**2 - 1.0) / 3.0
            valid_v = np.isfinite(E_ana_v)
            ax_v.plot(r_v, tau_v * Ev, ".", label=rf"Sim $\tau={tau_v:.2f}$")
            ax_v.plot(r_v[valid_v], tau_v * E_ana_v[valid_v], "-", alpha=0.7, label=rf"Ana $\tau={tau_v:.2f}$")

        ax_i.set_xlabel("r")
        ax_i.set_ylabel(r"$\tau T^{00}$")
        ax_i.set_title("Ideal")
        ax_i.grid(True)
        ax_i.legend(fontsize=8)

        ax_v.set_xlabel("r")
        ax_v.set_title("Viscous")
        ax_v.grid(True)
        ax_v.legend(fontsize=8)

        fig.suptitle("Ideal and Viscous: Sim vs Ana (Matched Times)")
        fig.tight_layout()
        out_overlay = os.path.join(OUTDIR, "gubser_ideal_viscous_matched.png")
        fig.savefig(out_overlay)
        print(f"Overlay plot saved to {out_overlay}")


        # Fractional difference plot: (E_visc - E_ideal) / E_ideal at matched times.
        fig_delta, ax_delta = plt.subplots(1, 1, figsize=(10, 6))
        xplot_v = np.abs(xs_v)
        order_v = np.argsort(xplot_v)

        for iv in v_indices:
            tau_v = float(tv[iv])
            ii = int(np.argmin(np.abs(ti - tau_v)))
            tau_i = float(ti[ii])

            Ev = fv["solution"][iv, mid_y_v, :, 0]
            Ei = fi["solution"][ii, mid_y_i, :, 0]

            # Interpolate ideal profile to viscous x-grid before differencing.
            Ei_on_v = np.interp(xs_v, xs_i, Ei)
            valid = np.abs(Ei_on_v) > 1e-14
            frac = np.full_like(Ev, np.nan, dtype=float)
            frac[valid] = (Ev[valid] - Ei_on_v[valid]) / Ei_on_v[valid]

            ax_delta.plot(
                xplot_v[order_v],
                frac[order_v],
                "-",
                label=rf"$\tau_v={tau_v:.2f}$ vs $\tau_i={tau_i:.2f}$",
            )

        ax_delta.axhline(0.0, color="k", linestyle="--", linewidth=1)
        ax_delta.set_xlabel(r"$r=|x|$ at mid-$y$")
        ax_delta.set_ylabel(r"$(E_{visc}-E_{ideal})/E_{ideal}$")
        ax_delta.set_title("Viscous vs Ideal Fractional Difference")
        ax_delta.grid(True)
        ax_delta.legend(fontsize=8)
        fig_delta.tight_layout()
        out_delta = os.path.join(OUTDIR, "gubser_ideal_viscous_fractional_diff.png")
        fig_delta.savefig(out_delta)
        print(f"Fractional-difference plot saved to {out_delta}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Viscous Gubser Flow Verification Test",
        usage="%(prog)s {run|plot} [options]",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the viscous Gubser flow verification test")
    run_parser.add_argument(
        "--run_command",
        required=True,
        type=str,
        default="./VischydroMain.exe",
        help="Command to run the VischydroMain executable",
    )
    run_parser.add_argument("--q", type=float, default=DEFAULT_Q, help="Gubser q [1/fm]")
    run_parser.add_argument("--eta_by_s", type=float, default=DEFAULT_ETA_BY_S, help="Shear viscosity to entropy density ratio")
    run_parser.add_argument("--hatT0", type=float, default=DEFAULT_HAT_T0, help="Dimensionless hat{T}_0 from paper")
    run_parser.add_argument("--f_star", type=float, default=DEFAULT_F_STAR, help="Conformal EOS constant f_*")
    run_parser.add_argument("--tau0", type=float, default=DEFAULT_TAU0, help="Initial proper time [fm/c]")
    run_parser.add_argument("--t_end", type=float, default=DEFAULT_T_END, help="Final proper time [fm/c]")
    run_parser.add_argument("--dt", type=float, default=0.075, help="Maximum time step")
    run_parser.add_argument("--print_frequency", type=int, default=1, help="Save every N solver steps")
    run_parser.add_argument("--nx", type=int, default=60, help="Grid points in x")
    run_parser.add_argument("--ny", type=int, default=60, help="Grid points in y")

    plot_parser = subparsers.add_parser("plot", help="Plot results from viscous Gubser flow test")
    plot_parser.add_argument("--q", type=float, default=DEFAULT_Q, help="Gubser q [1/fm]")
    plot_parser.add_argument("--eta_by_s", type=float, default=DEFAULT_ETA_BY_S, help="Shear viscosity to entropy density ratio")
    plot_parser.add_argument("--hatT0", type=float, default=DEFAULT_HAT_T0, help="Dimensionless hat{T}_0 from paper")
    plot_parser.add_argument("--f_star", type=float, default=DEFAULT_F_STAR, help="Conformal EOS constant f_*")
    plot_parser.add_argument("--max_curves", type=int, default=6, help="Maximum number of time slices to plot")
    plot_parser.add_argument(
        "--no_overlay_ideal",
        action="store_true",
        help="Disable auto-overlay with ideal run, even if ideal data files are found",
    )
    plot_parser.add_argument(
        "--ideal_run_name",
        type=str,
        default="gubser_test",
        help="Ideal run prefix used to load <run>_grid.h5 and <run>_grid_t.txt",
    )

    args = parser.parse_args()

    run_name = "gubser_viscous_test"

    if args.command == "run":
        run(
            run_name=run_name,
            q=args.q,
            hatT0=args.hatT0,
            tau0=args.tau0,
            eta_by_s=args.eta_by_s,
            f_star=args.f_star,
            t_end=args.t_end,
            dt=args.dt,
            print_frequency=args.print_frequency,
            nx=args.nx,
            ny=args.ny,
            exe_path=args.run_command,
        )
    elif args.command == "plot":
        analyze_results(
            run_name=run_name,
            q=args.q,
            hatT0=args.hatT0,
            eta_by_s=args.eta_by_s,
            f_star=args.f_star,
            max_curves=args.max_curves,
            overlay_ideal=(not args.no_overlay_ideal),
            ideal_run_name=args.ideal_run_name,
        )
