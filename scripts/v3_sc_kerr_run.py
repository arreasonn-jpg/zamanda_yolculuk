import numpy as np, time
from zt005.bssn_state_mp import BSSNStateMP
from zt005.sc_like_ic_kerr import sc_like_kerr_state
from zt005.bssn_rk4_mp import rk4_step_mp
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 20.0; M = 1.0; rc = 2.0; dt = 5e-4; steps = 50
h = L/(N-1)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R = np.sqrt(X**2 + Y**2 + Z**2)
mask = R > 3.0*rc

print(f"{'a':>5} {'H0':>11} {'Hf':>11} {'growth':>8} {'alpha_min':>10} {'bmax':>10}")
for a in [0.0, 0.1, 0.3, 0.5]:
    s = BSSNStateMP.from_base(sc_like_kerr_state(N, L, M=M, rc=rc, a=a))
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
    t0 = time.time()
    for _ in range(steps):
        s = rk4_step_mp(s, h=h, dt=dt)
    Hf = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
    print(f"{a:5.2f} {H0:11.4e} {Hf:11.4e} {Hf/H0:8.4f} "
          f"{s.alpha.min():10.4f} {np.max(np.abs(s.beta)):10.3e}")
