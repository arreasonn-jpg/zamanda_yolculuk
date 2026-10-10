import numpy as np, time
from zt005.bssn_ic import schwarzschild_isotropic_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 7; L = 20.0; M = 1.0; dt = 1e-3; n_steps = 500
h = L/(N-1)
s = schwarzschild_isotropic_state(N, L, M)
print(f"N={N} h={h:.2f} dt={dt:.1e} steps={n_steps}")
H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
    s.phi, s.gtilde, s.K, s.Atilde, h))))
print(f"H0 = {H0:.4e}")
norms = [H0]
t0 = time.time()
for step in range(1, n_steps+1):
    s = rk4_step_beta3(s, h=h, dt=dt)
    if step % 50 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        n = float(np.mean(np.abs(H)))
        norms.append(n)
        print(f"step={step:4d} H={n:.4e} growth={n/H0:.4f} "
              f"alpha_min={s.alpha.min():.4f} t={time.time()-t0:.0f}s")
print(f"final growth = {norms[-1]/norms[0]:.4e}")
print(f"all finite: {all(np.isfinite(norms))}")
