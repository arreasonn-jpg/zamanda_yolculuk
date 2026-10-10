import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_beta1 import rk4_step_beta1
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 7; L = 4.0; h = L/(N-1); dt = 1e-4
s = minkowski_state(N)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
s.phi = 1e-3 * np.exp(-(X**2 + Y**2 + Z**2) / 0.5)

norms = []
t0 = time.time()
for step in range(200):
    s = rk4_step_beta1(s, h=h, dt=dt)
    if step % 10 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        n = float(np.mean(np.abs(H)))
        norms.append(n)
        print(f"step={step:4d} H={n:.4e} alpha_min={s.alpha.min():.4f} Kmax={np.max(np.abs(s.K)):.2e} Atmax={np.max(np.abs(s.Atilde)):.2e}")
print(f"elapsed = {time.time()-t0:.1f}s")
print(f"H_first = {norms[0]:.4e}  H_last = {norms[-1]:.4e}")
print(f"total growth = {norms[-1]/norms[0]:.4e}")
print(f"all finite: {all(np.isfinite(norms))}")
