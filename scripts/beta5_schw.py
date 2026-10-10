import numpy as np, time
from zt005.bssn_ic import schwarzschild_isotropic_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 40.0; M = 1.0; h = L/(N-1); dt = 5e-3
s = schwarzschild_isotropic_state(N, L, M)
print(f"h={h:.3f} dt={dt:.2e} CFL={dt/h:.3f}")

H0 = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
print(f"H0 mean|.|={np.mean(np.abs(H0)):.4e} max={np.max(np.abs(H0)):.4e}")

norms = []
t0 = time.time()
for step in range(100):
    s = rk4_step_beta3(s, h=h, dt=dt)
    if step % 10 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        n = float(np.mean(np.abs(H)))
        norms.append(n)
        print(f"step={step:4d} H={n:.4e} alpha_min={s.alpha.min():.4f} t={time.time()-t0:.0f}s")
print(f"growth = {norms[-1]/norms[0]:.4e}")
print(f"all finite: {all(np.isfinite(norms))}")
