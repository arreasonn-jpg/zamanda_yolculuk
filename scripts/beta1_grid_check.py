import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_beta1 import rk4_step_beta1
from zt005.bssn_ricci_full import hamiltonian_constraint_full as hamiltonian_constraint_grid

def run(N, n_steps=5, dt=1e-4, A=1e-3, sigma=0.5):
    L = 4.0
    h = L / (N - 1)
    s = minkowski_state(N)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.phi = A * np.exp(-(X**2 + Y**2 + Z**2) / (2*sigma**2))
    norms = []
    t0 = time.time()
    for step in range(n_steps):
        s = rk4_step_beta1(s, h=h, dt=dt)
        H = hamiltonian_constraint_grid(s.phi, s.gtilde, s.K, s.Atilde, h)
        norms.append(float(np.mean(np.abs(H))))
    return {
        "N": N, "h": h, "steps": n_steps,
        "H0": norms[0], "Hf": norms[-1],
        "growth": norms[-1]/max(norms[0], 1e-30),
        "per_step": (norms[-1]/max(norms[0],1e-30))**(1.0/n_steps),
        "elapsed": time.time()-t0,
        "all_finite": all(np.isfinite(norms)),
    }

for N in [7, 9, 11]:
    r = run(N)
    print(f"N={r['N']:2d} h={r['h']:.3f} H0={r['H0']:.2e} Hf={r['Hf']:.2e} "
          f"growth={r['growth']:.2e} per_step={r['per_step']:.3f} t={r['elapsed']:.1f}s")
