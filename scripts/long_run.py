import numpy as np, time
from zt005.bssn_stability import gaussian_phi_grid
from zt005.bssn_rk4_grid import rk4_step_grid
from zt005.bssn import BSSN, hamiltonian_constraint

N = 11; L = 4.0; A = 1e-4; sigma = 0.5
h = L / (N - 1)
state = {
    "phi": gaussian_phi_grid(N, L, A, sigma),
    "gt":  np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy(),
    "K":   np.zeros((N,N,N)),
    "At":  np.zeros((N,N,N,3,3)),
}
norms = []
t0 = time.time()
for step in range(20):
    state = rk4_step_grid(state, h=h, dt=1e-4)
    H = 0.0
    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                b = BSSN(float(state["phi"][ix,iy,iz]), state["gt"][ix,iy,iz],
                         float(state["K"][ix,iy,iz]), state["At"][ix,iy,iz], np.zeros(3))
                H += abs(hamiltonian_constraint(b, 0.0))
    norms.append(H / N**3)
    if step % 5 == 0:
        print(f"step={step:3d} meanH={norms[-1]:.3e} maxphi={np.max(np.abs(state['phi'])):.3e}")
print(f"elapsed = {time.time()-t0:.1f}s")
print(f"final maxH = {norms[-1]:.3e}  maxphi = {np.max(np.abs(state['phi'])):.3e}")
print(f"H bounded: {all(np.isfinite(norms))}  last/first = {norms[-1]/max(norms[0],1e-30):.3e}")
