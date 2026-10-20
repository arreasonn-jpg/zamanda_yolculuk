"""v7.4: detailed study of QFT-compatible sub-regime (A < 1.87).

For A in [0.5, 1.87], measure:
- rho_min scaling with A
- NEC violation strength
- Whether PTUM's "discharge" is still meaningful at low amplitude
"""
import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex_ptum import rk4_step_complex_ptum
from zt005.matter.ptum_field import psi_equilibrium_smooth
from zt005.matter.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 100
params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

x = np.linspace(-L/2, L/2, N)
X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z2**2

A_list = [0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
print(f"Sub-QFT regime: A in [0.5, 1.8]")
print(f"{'A':>5} {'rho_min':>12} {'ratio':>10} {'n_neg':>6} "
      f"{'Phimax':>8} {'psimax':>8} {'viol':>6}")

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
    n_neg = int(np.sum(rho_arr < 0))
    Phimax = float(np.max(np.abs(Phi)))
    psimax = float(np.max(np.abs(psi)))
    print(f"{A:>5.2f} {rho_min:>12.4e} {r['ratio']:>10.4f} {n_neg:>6d} "
          f"{Phimax:>8.4f} {psimax:>8.4f} {str(r['violates']):>6} "
          f"t={time.time()-t0:.0f}s")
    results.append((A, rho_min, r['ratio'], n_neg, Phimax, psimax))

np.savez("results/validated/v7_4_sub_qft.npz",
         A=np.array([r[0] for r in results]),
         rho_min=np.array([r[1] for r in results]),
         ratio=np.array([r[2] for r in results]),
         n_neg=np.array([r[3] for r in results]))
print("Saved: results/validated/v7_4_sub_qft.npz")
