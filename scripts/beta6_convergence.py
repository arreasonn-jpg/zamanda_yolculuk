import numpy as np, time
from zt005.bssn_ic import schwarzschild_isotropic_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

L = 40.0; M = 1.0; dt = 5e-3
print(f"L={L} M={M} dt={dt}")
print(f"{'N':>3} {'h':>7} {'H0_mean':>11} {'H_10_mean':>11} {'ratio':>8} {'order':>6}")
results = []
for N in [7, 9, 11, 13]:
    h = L / (N - 1)
    s = schwarzschild_isotropic_state(N, L, M)
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    for _ in range(10):
        s = rk4_step_beta3(s, h=h, dt=dt)
    H10 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h))))
    results.append((N, h, H0, H10))
    print(f"{N:3d} {h:7.3f} {H0:11.4e} {H10:11.4e}", end="")
    if len(results) > 1:
        N_prev, h_prev, H0_prev, _ = results[-2]
        ratio = H0 / H0_prev
        order = np.log(ratio) / np.log(h / h_prev)
        print(f" {ratio:8.3f} {order:6.2f}")
    else:
        print()
