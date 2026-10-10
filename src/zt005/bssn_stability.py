"""3D stability test: small Gaussian perturbation on Minkowski, evolve
with bssn_rk4, verify constraint norm stays bounded."""
from __future__ import annotations
import numpy as np
from .bssn import BSSN, hamiltonian_constraint
from .bssn_rk4 import evolve_rk4


def gaussian_phi_grid(N, L, A, sigma):
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2
    return A * np.exp(-R2 / (2.0 * sigma**2))


def run_stability(N=11, L=4.0, A=1e-4, sigma=0.5, dt=1e-3, n_steps=20):
    phi0 = gaussian_phi_grid(N, L, A, sigma)
    gt0 = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    K0 = np.zeros((N,N,N))
    At0 = np.zeros((N,N,N,3,3))
    h = L / (N - 1)
    out = evolve_rk4(phi0, gt0, K0, At0, h=h, dt=dt, n_steps=n_steps)
    # constraint norm over final state
    H = []
    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                b = BSSN(float(out["phi"][ix,iy,iz]), out["gtilde"][ix,iy,iz],
                         float(out["K"][ix,iy,iz]), out["Atilde"][ix,iy,iz],
                         np.zeros(3))
                H.append(abs(hamiltonian_constraint(b, 0.0)))
    return {"phi": out["phi"], "constraint_norm": np.array(H),
            "max_phi": float(np.max(np.abs(out["phi"]))),
            "grid_h": h, "dt": dt, "n_steps": n_steps}


__all__ = ["gaussian_phi_grid", "run_stability"]
