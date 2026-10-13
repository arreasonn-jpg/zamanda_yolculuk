"""v7.3-alpha2: fine threshold scan + paper figure."""
import numpy as np, time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex_ptum import rk4_step_complex_ptum
from zt005.ptum_field import psi_equilibrium_smooth
from zt005.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 100
params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

x = np.linspace(-L/2, L/2, N)
X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z2**2

A_list = [1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.5, 3.0, 3.5]
results = []

for A in A_list:
    t0 = time.time()
    s = minkowski_state(N)
    Z = np.zeros((N,N,N))
    Phi = A*np.exp(-R2/8.0)
    psi = psi_equilibrium_smooth(Phi)
    PR = 1e-4*np.exp(-R2/2.0)
    PI = 0.5e-4*np.exp(-R2/2.0)

    t_hist = []; rho_hist = []
    for step in range(steps):
        s, PR, PI, Z, Z, Phi, Z, psi, Z = rk4_step_complex_ptum(
            s, PR, PI, Z, Z, Phi, Z, psi, Z, h, dt, params)
        if (step + 1) % 5 == 0:
            src = ptum_sources(Phi, Z, psi, Z, s.gtilde, s.phi, h, **ptum_kwargs)
            t_hist.append((step+1)*dt)
            rho_hist.append(float(src["rho"].min()))
    t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
    rho_min = float(np.min(rho_arr))
    r = check_ford_roman(rho_arr, t_arr, tau=1.0, N_field_components=2)
    results.append((A, rho_min, r["ratio"], r["violates"]))
    print(f"A={A:.2f} rho_min={rho_min:.4e} ratio={r['ratio']:.4f} "
          f"viol={r['violates']} t={time.time()-t0:.0f}s")

A_arr = np.array([r[0] for r in results])
rho_arr = np.array([r[1] for r in results])
ratio_arr = np.array([r[2] for r in results])

np.savez("results/validated/v7_3_qft_scan.npz",
         A=A_arr, rho_min=rho_arr, ratio=ratio_arr)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].plot(np.abs(rho_arr), ratio_arr, 'bo-', markersize=6)
axes[0].axhline(1.0, color='r', ls='--', label='QFT threshold (ratio=1)')
A_fit = np.polyfit(np.abs(rho_arr), ratio_arr, 1)
x_fit = np.linspace(0, np.max(np.abs(rho_arr)), 50)
axes[0].plot(x_fit, A_fit[0]*x_fit + A_fit[1], 'g--',
              label=f'fit: {A_fit[0]:.3f}*|rho| + {A_fit[1]:.3f}')
axes[0].set_xlabel(r'$|\rho_{\min}|$')
axes[0].set_ylabel('Ford-Roman ratio (tau=1)')
axes[0].set_title('(a) QFT violation vs NEC magnitude')
axes[0].legend(); axes[0].grid(True, alpha=0.3)

axes[1].plot(A_arr, ratio_arr, 'mo-', markersize=6)
axes[1].axhline(1.0, color='r', ls='--')
axes[1].set_xlabel(r'$A_{\rm fuel}$')
axes[1].set_ylabel('Ford-Roman ratio (tau=1)')
axes[1].set_title('(b) QFT violation vs fuel amplitude')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("paper/figures/fig8_qft_threshold.png", dpi=100)
print("Saved: paper/figures/fig8_qft_threshold.png")
print(f"Linear fit: ratio = {A_fit[0]:.4f} * |rho_min| + {A_fit[1]:.4f}")
