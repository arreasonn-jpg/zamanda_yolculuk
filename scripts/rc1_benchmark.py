import numpy as np, time
from zt005.bssn_ic import schwarzschild_isotropic_state
from zt005.bssn_ic_kerr import kerr_slow_rotation_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 7; L = 20.0; M = 1.0; dt = 1e-3; steps = 100
h = L/(N-1)
print(f"{'case':>12} {'H0':>11} {'Hf':>11} {'growth':>8} {'alpha_min':>10} {'bmax':>10}")
cases = [("Schwarz", 0.0, schwarzschild_isotropic_state)]
for a in [0.0, 0.2, 0.4, 0.5]:
    cases.append((f"Kerr a={a}", a, None))

for name, a, ic_fn in cases:
    if ic_fn is not None:
        s = ic_fn(N, L, M)
    else:
        s = kerr_slow_rotation_state(N, L, M=M, a=a)
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    for _ in range(steps):
        s = rk4_step_beta3(s, h=h, dt=dt)
    Hf = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    print(f"{name:>12} {H0:11.4e} {Hf:11.4e} {Hf/H0:8.4f} "
          f"{s.alpha.min():10.4f} {np.max(np.abs(s.beta)):10.3e}")
