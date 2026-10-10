import numpy as np, time
from zt005.bssn_ic_kerr_full import kerr_schild_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 40.0; M = 1.0; dt = 5e-4; steps = 50
h = L/(N-1)
print(f"N={N} L={L} h={h:.2f} dt={dt} steps={steps}")
print(f"{'a':>5} {'H0':>11} {'Hf':>11} {'growth':>8} {'alpha_min':>10} {'t':>5}")
for a in [0.0, 0.3, 0.6, 0.9, 0.99]:
    s = kerr_schild_state(N, L, M=M, a=a)
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    t0 = time.time()
    for _ in range(steps):
        s = rk4_step_beta3(s, h=h, dt=dt)
    Hf = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    print(f"{a:5.2f} {H0:11.4e} {Hf:11.4e} {Hf/H0:8.4f} "
          f"{s.alpha.min():10.4f} {time.time()-t0:5.0f}")
