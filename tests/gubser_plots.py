import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# --- Parameters for Au-Au Collision ---
q = 1 / 4.3           # Inverse length scale (1/fm)
T0_hat = 6.0          # de Sitter temperature scale (dimensionless)
r_max = 12.0          # Transverse radius (fm)
points = 500          # Number of radial points
tau_vals = [0.6, 1.5, 3.0, 5.0, 8.0]  # Proper times in fm/c
T_c = 155             # Phase transition temperature in MeV

# Scaling factor to convert T_hat/tau to MeV
# Usually related to hbar*c if using natural units, 
# here we treat T0_hat such that peak T at tau=0.6 is ~400-500 MeV.
unit_conv = 197.3  # Converting 1/fm to MeV

def get_rho(tau, r, q):
    """Maps Minkowski (tau, r) to de Sitter rho."""
    return np.arcsinh((1 - q**2 * tau**2 + q**2 * r**2) / (2 * q * tau))

def ideal_gubser_ode(T_hat, rho):
    """The governing ODE: d(T_hat)/drho = -2/3 * T_hat * tanh(rho)"""
    return -(2/3) * T_hat * np.tanh(rho)

def transform_to_minkowski(T_hat_val, rho_val, tau):
    """T(tau, r) = T_hat(rho) / (tau * cosh(rho))"""
    return (T_hat_val * unit_conv) / (tau * np.cosh(rho_val))

# --- STEP 1: Numerical Master Solution ---
rho_limit = 10.0
rho_pos = np.linspace(0, rho_limit, 4000)
sol_pos = odeint(ideal_gubser_ode, T0_hat, rho_pos)
rho_neg = np.linspace(0, -rho_limit, 4000)
sol_neg = odeint(ideal_gubser_ode, T0_hat, rho_neg)
rho_full = np.concatenate([rho_neg[::-1], rho_pos[1:]])
T_hat_full = np.concatenate([sol_neg[::-1], sol_pos[1:]]).flatten()

# --- STEP 2: Plotting ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
colors = plt.cm.plasma(np.linspace(0, 0.85, len(tau_vals)))

for i, tau in enumerate(tau_vals):
    r_vals = np.linspace(0, r_max, points)
    rho_vals = get_rho(tau, r_vals, q)
    
    # 1. Temperature (MeV)
    T_hat_num = np.interp(rho_vals, rho_full, T_hat_full)
    T_num = transform_to_minkowski(T_hat_num, rho_vals, tau)
    T_analyt = transform_to_minkowski(T0_hat / (np.cosh(rho_vals)**(2/3)), rho_vals, tau)
    
    # 2. Energy Density (epsilon ~ T^4)
    # Factor 1e-8 is a scaling for visibility; 
    # in physical units e = (g_pi^2/30) * T^4
    e_num = T_num**4 / 1e8 
    e_analyt = T_analyt**4 / 1e8
    
    # 3. Transverse Velocity
    v_perp = (2 * q**2 * tau * r_vals) / (1 + q**2 * tau**2 + q**2 * r_vals**2)
    
    # Subplot 1: Temperature (Linear plot)
    axes[0].plot(r_vals, T_analyt, color=colors[i], lw=2, label=rf'$\tau={tau}$ fm/c')
    axes[0].plot(r_vals, T_num, 'k--', lw=0.8, alpha=0.5)
    
    # Subplot 2: Energy Density (Log plot)
    axes[1].semilogy(r_vals, e_analyt, color=colors[i], lw=2)
    axes[1].semilogy(r_vals, e_num, 'k--', lw=0.8, alpha=0.5)
    
    # Subplot 3: Velocity (Linear plot)
    axes[2].plot(r_vals, v_perp, color=colors[i], lw=2)

# Styling
axes[0].axhline(T_c, color='gray', linestyle=':', label='$T_c = 155$ MeV')
axes[0].set_title('Temperature $T$ [MeV]', fontsize=14)
axes[0].set_ylim(0, 600)

axes[1].set_title('Energy Density $\epsilon$ [arb. units]', fontsize=14)
axes[1].set_ylim(1e-2, 2e3)

axes[2].set_title('Transverse Velocity $v_\perp$', fontsize=14)
axes[2].set_ylim(0, 1.0)

for ax in axes:
    ax.set_xlabel(r'$r$ [fm]', fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.set_xlim(0, r_max)

axes[0].legend(loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig('gubser_ideal_plots.png', dpi=300)
plt.show()