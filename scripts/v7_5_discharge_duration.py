"""v7.5-alpha1: discharge duration in the QFT-compatible window.

For A in {1.4, 1.6, 1.8}, run 500 steps and measure:
- Phi decay time (Phi < 1.0)
- rho_min decay (rho > -0.1)
- Ford-Roman evolution with time window
"""
import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex_ptum import rk4_step_complex_ptum
from zt005.ptum_field import psi_equilibrium_smooth
from zt005.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 500
params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

x = np.linspace(-L/2, L/2, N)
X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z2**2

print(f"QFT-OK discharge duration: N={N} dt={dt} steps={steps}")

for A in [1.4, 1.6, 1.8]:
    t0 = time.time()
    s = minkowski_state(N)
    Z = np.zeros((N,N,N))
    Phi = A*np.exp(-R2/8.0)
    psi = psi_equilibrium_smooth(Phi)
    PR = 1e-4*np.exp(-R2/2.0)
    PI = 0.5e-4*np.exp(-R2/2.0)

    t_hist = []; rho_hist = []; phimax_hist = []
    for step in range(steps):
        s, PR, PI, Z, Z, Phi, Z, psi, Z = rk4_step_complex_ptum(
            s, PR, PI, Z, Z, Phi, Z, psi, Z, h, dt, params)
        if (step + 1) % 10 == 0:
            src = ptum_sources(Phi, Z, psi, Z, s.gtilde, s.phi, h, **ptum_kwargs)
            t_hist.append((step+1)*dt)
            rho_hist.append(float(src["rho"].min()))
            phimax_hist.append(float(np.max(np.abs(Phi))))

    t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
    phimax_arr = np.array(phimax_hist)

    # Discharge duration metrics
    phi_below = np.where(phimax_arr < 1.0)[0]
    phi_decay_t = t_hist[phi_below[0]] if len(phi_below) > 0 else None
    rho_above = np.where(rho_arr > -0.1)[0]
    rho_decay_t = t_hist[rho_above[0]] if len(rho_above) > 0 else None

    print(f"\nA={A}: t_wall={time.time()-t0:.0f}s")
    print(f"  rho_min: t=0 {rho_arr[0]:.4e}, t=end {rho_arr[-1]:.4e}")
    print(f"  Phimax:  t=0 {phimax_arr[0]:.4f}, t=end {phimax_arr[-1]:.4f}")
    if phi_decay_t is not None:
        print(f"  Phi < 1 at t = {phi_decay_t:.3f}")
    else:
        print(f"  Phi < 1: not reached in {steps*dt:.2f}")
    if rho_decay_t is not None:
        print(f"  rho > -0.1 at t = {rho_decay_t:.3f}")
    else:
        print(f"  rho > -0.1: not reached in {steps*dt:.2f}")
