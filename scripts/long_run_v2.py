import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_v2 import rk4_step_v2
from zt005.bssn import BSSN, hamiltonian_constraint

N = 7; L = 4.0
h = L / (N - 1)
s = minkowski_state(N)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
s.phi = 1e-5 * np.exp(-(X**2 + Y**2 + Z**2) / 0.5)

norms = []
t0 = time.time()
for step in range(40):
    s = rk4_step_v2(s, h=h, dt=2e-5)
    H = 0.0
    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                b = BSSN(float(s.phi[ix,iy,iz]), s.gtilde[ix,iy,iz],
                         float(s.K[ix,iy,iz]), s.Atilde[ix,iy,iz], np.zeros(3))
                H += abs(hamiltonian_constraint(b, 0.0))
    norms.append(H / N**3)
    if step % 5 == 0:
        print(f"step={step:3d} meanH={norms[-1]:.3e} alpha_min={s.alpha.min():.3f} Gmax={np.max(np.abs(s.Gamma)):.2e}")
print(f"elapsed = {time.time()-t0:.1f}s")
print(f"final meanH = {norms[-1]:.3e}")
print(f"growth last/first = {norms[-1]/max(norms[0],1e-30):.3e}")
print(f"all finite: {all(np.isfinite(norms))}")
