import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 7; L = 4.0; h = L/(N-1); dt = 5e-5
s = minkowski_state(N)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
s.phi = 1e-3 * np.exp(-(X**2 + Y**2 + Z**2) / 0.5)

norms = []
t0 = time.time()
for step in range(100):
    s = rk4_step_beta3(s, h=h, dt=dt)
    if step % 10 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        n = float(np.mean(np.abs(H)))
        norms.append(n)
        print(f"step={step:4d} H={n:.4e} alpha_min={s.alpha.min():.4f} "
              f"bmax={np.max(np.abs(s.beta)):.2e}")
print(f"elapsed = {time.time()-t0:.1f}s")
print(f"growth = {norms[-1]/norms[0]:.4e}")
print(f"all finite: {all(np.isfinite(norms))}")
