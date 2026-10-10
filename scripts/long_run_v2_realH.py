import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_v2 import rk4_step_v2
from zt005.bssn_constraints import hamiltonian_constraint_grid

N = 11; L = 4.0
h = L / (N - 1)
s = minkowski_state(N)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
s.phi = 1e-3 * np.exp(-(X**2 + Y**2 + Z**2) / 0.5)

norms = []
t0 = time.time()
for step in range(10):
    s = rk4_step_v2(s, h=h, dt=1e-5)
    H = hamiltonian_constraint_grid(s.phi, s.gtilde, s.K, s.Atilde, h)
    norms.append(float(np.mean(np.abs(H))))
    if step % 5 == 0:
        print(f"step={step:3d} meanH={norms[-1]:.3e} alpha_min={s.alpha.min():.4f} Kmax={np.max(np.abs(s.K)):.2e}")
print(f"elapsed = {time.time()-t0:.1f}s")
print(f"growth = {norms[-1]/max(norms[0],1e-30):.3e}")
print(f"all finite: {all(np.isfinite(norms))}")
