"""Grid-based constraint residuals with true Ricci scalar (v2.0-alpha6)."""
from __future__ import annotations
import numpy as np
from .bssn_grid import christoffel_from_grid, dG_from_grid
from .bssn_rhs import ricci_tensor


def hamiltonian_constraint_grid(phi, gt, K, At, h, rho_grid=None):
    """H = R + K^2 - K_ij K^ij - 16 pi G rho / c^2 (true Ricci scalar)."""
    C2 = (2.99792458e8) ** 2
    G = 6.67430e-11
    Gc, _ = christoffel_from_grid(gt, h)
    dG = dG_from_grid(Gc, h)
    Nx, Ny, Nz = phi.shape
    H = np.zeros_like(phi)
    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                R = ricci_tensor(Gc[ix,iy,iz], dG[ix,iy,iz])
                Rscalar = float(np.trace(R))
                A2 = float(np.sum(At[ix,iy,iz] ** 2))
                rho = float(rho_grid[ix,iy,iz]) if rho_grid is not None else 0.0
                H[ix,iy,iz] = (Rscalar + K[ix,iy,iz]**2 - A2
                               - 16.0 * np.pi * G * rho / C2)
    return H


__all__ = ["hamiltonian_constraint_grid"]
