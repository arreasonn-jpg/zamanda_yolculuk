"""RK4 driver for the BSSN skeleton (uses bssn_rhs + grid FD)."""
from __future__ import annotations
import numpy as np
from .bssn import BSSN
from .bssn_rhs import bssn_rhs
from .bssn_grid import dgtilde_from_grid, christoffel_from_grid, dG_from_grid


def _rhs_field(phi, gt, K, At, h, rho_fn, S_fn, eta):
    Nx, Ny, Nz = phi.shape
    dg = dgtilde_from_grid(gt, h)
    G, _ = christoffel_from_grid(gt, h)
    dG = dG_from_grid(G, h)
    dphi = np.zeros_like(phi); dgt = np.zeros_like(gt)
    dK = np.zeros_like(K); dAt = np.zeros_like(At)
    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                rho = rho_fn(ix,iy,iz) if rho_fn else 0.0
                S = S_fn(ix,iy,iz) if S_fn else np.zeros(3)
                b = BSSN(float(phi[ix,iy,iz]), gt[ix,iy,iz],
                         float(K[ix,iy,iz]), At[ix,iy,iz], np.zeros(3))
                r = bssn_rhs(b, rho, S, dg[ix,iy,iz], dG[ix,iy,iz], eta=eta)
                dphi[ix,iy,iz] = r["dphi"]
                dgt[ix,iy,iz]  = r["dgt"]
                dK[ix,iy,iz]   = r["dK"]
                dAt[ix,iy,iz]  = r["dAt"]
    return dphi, dgt, dK, dAt


def evolve_rk4(phi0, gt0, K0, At0, h, dt, n_steps, rho_fn=None,
               S_fn=None, eta=0.1):
    phi, gt, K, At = phi0.copy(), gt0.copy(), K0.copy(), At0.copy()
    norms = []
    for _ in range(n_steps):
        k1 = _rhs_field(phi, gt, K, At, h, rho_fn, S_fn, eta)
        k2 = _rhs_field(phi+0.5*dt*k1[0], gt+0.5*dt*k1[1],
                        K+0.5*dt*k1[2], At+0.5*dt*k1[3], h, rho_fn, S_fn, eta)
        k3 = _rhs_field(phi+0.5*dt*k2[0], gt+0.5*dt*k2[1],
                        K+0.5*dt*k2[2], At+0.5*dt*k2[3], h, rho_fn, S_fn, eta)
        k4 = _rhs_field(phi+dt*k3[0], gt+dt*k3[1],
                        K+dt*k3[2], At+dt*k3[3], h, rho_fn, S_fn, eta)
        phi = phi + (dt/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        gt  = gt  + (dt/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        K   = K   + (dt/6.0)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
        At  = At  + (dt/6.0)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
        norms.append(float(np.max(np.abs(phi)) + np.max(np.abs(K))))
    return {"phi": phi, "gtilde": gt, "K": K, "Atilde": At,
            "constraint_norm": np.array(norms)}


__all__ = ["evolve_rk4"]
