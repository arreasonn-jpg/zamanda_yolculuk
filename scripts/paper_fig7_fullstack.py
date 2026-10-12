"""Generate Figure 7: full-stack evolution diagnostics (for paper)."""
import numpy as np, time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.kerr_sc_ptum_ic import kerr_sc_ptum_state
from zt005.bssn_rk4_sc_ptum import rk4_step_sc_ptum
from zt005.matter_ptum import ptum_sources

N = 9; L = 20.0; h = L/(N-1); dt = 1e-4; steps = 100
a = 0.9; A_fuel = 1.5
print(f"Generating full-stack figure: a={a} A_fuel={A_fuel} steps={steps}")

pk = kerr_sc_ptum_state(N=N, L=L, M=1.0, a=a, rc=2.0,
                         A_fuel=A_fuel, sigma_fuel=4.0)
s = pk["state"]
phi_sc = pk["phi_sc"]; Pi_sc = pk["Pi_sc"]
Phi = pk["Phi"]; Pi_Phi = pk["Pi_Phi"]
psi = pk["psi"]; Pi_psi = pk["Pi_psi"]
V_sc = pk["V_sc"]
params = pk["V_ptum_params"]

steps_arr = []
rho_arr = []
Kmax_arr = []
alpha_arr = []

for step in range(steps):
    s, phi_sc, Pi_sc, Phi, Pi_Phi, psi, Pi_psi = rk4_step_sc_ptum(
        s, phi_sc, Pi_sc, Phi, Pi_Phi, psi, Pi_psi,
        h, dt, V_sc, params)
    if (step + 1) % 5 == 0:
        ptum_kwargs = dict(params)
        if "alpha_p" in ptum_kwargs:
            ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")
        src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi,
                           s.gtilde, s.phi, h, **ptum_kwargs)
        steps_arr.append(step + 1)
        rho_arr.append(float(src["rho"].min()))
        Kmax_arr.append(float(np.max(np.abs(s.K))))
        alpha_arr.append(float(s.alpha.min()))

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

axes[0].plot(steps_arr, rho_arr, 'b-o', markersize=4)
axes[0].axhline(0, color='k', lw=0.5, ls='--')
axes[0].set_xlabel('Step'); axes[0].set_ylabel(r'$\rho_{\min}$')
axes[0].set_title('(a) Negative energy density')
axes[0].grid(True, alpha=0.3)

axes[1].plot(steps_arr, Kmax_arr, 'r-s', markersize=4)
axes[1].set_xlabel('Step'); axes[1].set_ylabel(r'$K_{\max}$')
axes[1].set_title('(b) Extrinsic curvature')
axes[1].grid(True, alpha=0.3)

axes[2].plot(steps_arr, alpha_arr, 'g-^', markersize=4)
axes[2].set_xlabel('Step'); axes[2].set_ylabel(r'$\alpha_{\min}$')
axes[2].set_title('(c) Lapse minimum')
axes[2].grid(True, alpha=0.3)

plt.suptitle(f'Full stack: Kerr+SC+condensate+PTUM (a={a}, A_fuel={A_fuel})',
             fontsize=11)
plt.tight_layout()
plt.savefig("paper/figures/fig7_fullstack.png", dpi=100)
print("Saved: paper/figures/fig7_fullstack.png")
print(f"rho_min: {rho_arr[0]:.4e} -> {rho_arr[-1]:.4e}")
print(f"Kmax:    {Kmax_arr[0]:.4e} -> {Kmax_arr[-1]:.4e}")
print(f"alpha:   {alpha_arr[0]:.4f} -> {alpha_arr[-1]:.4f}")
