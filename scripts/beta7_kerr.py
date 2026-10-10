import numpy as np, time
from zt005.bssn_ic_kerr import kerr_slow_rotation_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 20.0; M = 1.0; dt = 5e-4
print(f"N={N} L={L} dt={dt}")
print(f"{'a':>5} {'H0':>11} {'H50':>11} {'growth':>8} {'alpha_min':>10} {'bmax':>10}")
for a in [0.0, 0.1, 0.3, 0.5]:
    s = kerr_slow_rotation_state(N, L, M=M, a=a)
    h = L/(N-1)
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    t0 = time.time()
    for _ in range(50):
        s = rk4_step_beta3(s, h=h, dt=dt)
    H50 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    print(f"{a:5.2f} {H0:11.4e} {H50:11.4e} {H50/H0:8.3f} "
          f"{s.alpha.min():10.4f} {np.max(np.abs(s.beta)):10.3e} t={time.time()-t0:.0f}s")
