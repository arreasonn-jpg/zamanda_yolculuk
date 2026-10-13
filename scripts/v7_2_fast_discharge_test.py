"""v7.2: fast-discharge variant. Higher lambda -> faster Phi decay.

NOTE: lambda=1.0 is NOT the PTUM fiducial (lambda=0.05). This is
an explicitly labeled "fast-discharge" variant to test whether
transient NEC violation can be QFT-compatible.
"""
import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex_ptum import rk4_step_complex_ptum
from zt005.ptum_field import fuel_profile, psi_equilibrium_smooth
from zt005.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman, check_ford_roman_massive

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 500
A_fuel = 3.0
print(f"FAST DISCHARGE: A={A_fuel} lambda=1.0 (not fiducial)")

# Fast-discharge parameters
params = {"V0": 1.0, "lam": 1.0, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

s = minkowski_state(N)
Z = np.zeros((N,N,N))
Phi = fuel_profile(N, L, Phi0=A_fuel, sigma=4.0)
psi = psi_equilibrium_smooth(Phi)
PR = np.zeros((N,N,N)); PI = np.zeros((N,N,N))

t_hist = []; rho_hist = []; phimax_hist = []
t0 = time.time()
for step in range(steps):
    s, PR, PI, Z, Z, Phi, Z, psi, Z = rk4_step_complex_ptum(
        s, PR, PI, Z, Z, Phi, Z, psi, Z, h, dt, params)
    if (step + 1) % 25 == 0:
        src = ptum_sources(Phi, Z, psi, Z, s.gtilde, s.phi, h, **ptum_kwargs)
        t_hist.append((step+1)*dt)
        rho_hist.append(float(src["rho"].min()))
        phimax_hist.append(float(np.max(np.abs(Phi))))
        if (step + 1) % 100 == 0:
            print(f"step={step+1:4d} t={(step+1)*dt:.3f} "
                  f"rho_min={rho_hist[-1]:.4e} Phimax={phimax_hist[-1]:.4f} "
                  f"psimax={np.max(np.abs(psi)):.4f} t_wall={time.time()-t0:.0f}s")

# Save
np.savez("results/validated/v7_2_fast_discharge.npz",
         t=np.array(t_hist), rho=np.array(rho_hist),
         phimax=np.array(phimax_hist))

t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
print("\n=== Ford-Roman (massless, N=2) ===")
print(f"{'tau':>8} {'value':>14} {'bound':>14} {'ratio':>10} {'viol':>6}")
for tau in [0.01, 0.05, 0.1, 0.5, 1.0]:
    r = check_ford_roman(rho_arr, t_arr, tau=tau, N_field_components=2)
    print(f"{tau:>8.3f} {r['value']:>14.4e} {r['bound']:>14.4e} "
          f"{r['ratio']:>10.4f} {str(r['violates']):>6}")

# When does Phi decay below 1?
below = np.array(phimax_hist) < 1.0
if below.any():
    idx = np.argmax(below)
    print(f"\nPhi < 1 at t = {t_hist[idx]:.3f}")
below2 = rho_arr > -1.0
if below2.any():
    idx2 = np.argmax(below2)
    print(f"rho > -1 at t = {t_hist[idx2]:.3f}")
