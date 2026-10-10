import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full

# Domain kucuk: dalga sinirlara ulassin
N = 7; L = 0.5; h = L/(N-1); dt = 1e-3
print(f"h={h:.4f} dt={dt:.1e} CFL={dt/h:.3f} (must be <1)")
s = minkowski_state(N)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
s.phi = 1e-3 * np.exp(-(X**2 + Y**2 + Z**2) / 0.01)

norms = []
t0 = time.time()
for step in range(300):
    s = rk4_step_beta3(s, h=h, dt=dt)
    if step % 30 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        n = float(np.mean(np.abs(H)))
        norms.append(n)
        # boundary and center values
        bc_edge = abs(s.phi[0,0,0]) + abs(s.phi[-1,-1,-1])
        center = abs(s.phi[N//2,N//2,N//2])
        print(f"step={step:4d} H={n:.3e} edge={bc_edge:.2e} center={center:.2e} t={time.time()-t0:.0f}s")
print(f"growth = {norms[-1]/norms[0]:.4e}")
print(f"all finite: {all(np.isfinite(norms))}")
