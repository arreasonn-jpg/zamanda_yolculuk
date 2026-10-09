"""Time loop for BSSN skeleton on a uniform grid (Euler, no RK)."""
from __future__ import annotations
import numpy as np
from .bssn import BSSN, euler_step, hamiltonian_constraint
from .bssn_rhs import bssn_rhs
from .bssn_grid import dgtilde_from_grid, christoffel_from_grid, dG_from_grid


def evolve_bssn(phi0, gtilde0, K0, Atilde0, h, dt, n_steps,
                rho_fn=None, S_fn=None, eta=0.1):
    """Pointwise Euler time loop over an (Nx,Ny,Nz) grid.
    All state arrays share the same grid shape; derivatives via FD."""
    Nx, Ny, Nz = phi0.shape
    phi = phi0.copy(); gt = gtilde0.copy(); K = K0.copy(); At = Atilde0.copy()
    constraints = []
    for step in range(n_steps):
        dg = dgtilde_from_grid(gt, h)
        G, _ = christoffel_from_grid(gt, h)
        dG = dG_from_grid(G, h)
        for ix in range(Nx):
            for iy in range(Ny):
                for iz in range(Nz):
                    rho = rho_fn(ix, iy, iz) if rho_fn else 0.0
                    S = S_fn(ix, iy, iz) if S_fn else np.zeros(3)
                    b = BSSN(phi=float(phi[ix,iy,iz]),
                             gtilde=gt[ix,iy,iz],
                             K=float(K[ix,iy,iz]),
                             Atilde=At[ix,iy,iz],
                             Gamma=np.zeros(3))
                    rhs = bssn_rhs(b, rho, S, dg[ix,iy,iz], dG[ix,iy,iz], eta=eta)
                    b2 = euler_step(b, rhs, dt)
                    phi[ix,iy,iz] = b2.phi
                    gt[ix,iy,iz] = b2.gtilde
                    K[ix,iy,iz] = b2.K
                    At[ix,iy,iz] = b2.Atilde
        # monitor mean |H|
        Hs = []
        for ix in range(Nx):
            for iy in range(Ny):
                for iz in range(Nz):
                    b = BSSN(float(phi[ix,iy,iz]), gt[ix,iy,iz],
                             float(K[ix,iy,iz]), At[ix,iy,iz], np.zeros(3))
                    rho = rho_fn(ix,iy,iz) if rho_fn else 0.0
                    Hs.append(abs(hamiltonian_constraint(b, rho)))
        constraints.append(float(np.mean(Hs)))
    return {"phi": phi, "gtilde": gt, "K": K, "Atilde": At,
            "constraint_norm": np.array(constraints)}


__all__ = ["evolve_bssn"]
