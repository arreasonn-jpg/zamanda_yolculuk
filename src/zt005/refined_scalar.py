"""Two-grid scalar field evolution (v6.0-alpha4).

Coarse grid: N0, h0, domain [-L0/2, L0/2]^3
Fine grid:   N1, h1 = h0/2, centered at origin, domain [-L1/2, L1/2]^3
             with L1 = L0/2 (typically)

Coupling per step:
1. Evolve coarse RK4 (Sommerfeld BC on outer edges)
2. Prolongate coarse -> fine boundary values
3. Evolve fine RK4 (interior only)
4. Restrict fine -> coarse interior
"""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central
from .amr.refinement import prolongate, restrict


def _lap(f, h):
    grad = np.stack(d_central(f, h), axis=-1)
    L = np.zeros_like(f)
    for ax in range(3):
        dd = d_central(grad[..., ax], h)
        L += dd[ax]
    return L


def _scalar_rhs(phi, Pi, h, m2=0.0):
    """d_t phi = Pi; d_t Pi = Lap(phi) - m2 phi."""
    return Pi, _lap(phi, h) - m2 * phi


def _rk4_step(phi, Pi, h, dt, m2=0.0):
    d1phi, d1Pi = _scalar_rhs(phi, Pi, h, m2)
    p2 = phi + 0.5*dt*d1phi; P2 = Pi + 0.5*dt*d1Pi
    d2phi, d2Pi = _scalar_rhs(p2, P2, h, m2)
    p3 = phi + 0.5*dt*d2phi; P3 = Pi + 0.5*dt*d2Pi
    d3phi, d3Pi = _scalar_rhs(p3, P3, h, m2)
    p4 = phi + dt*d3phi; P4 = Pi + dt*d3Pi
    d4phi, d4Pi = _scalar_rhs(p4, P4, h, m2)
    new_phi = phi + (dt/6)*(d1phi + 2*d2phi + 2*d3phi + d4phi)
    new_Pi  = Pi  + (dt/6)*(d1Pi  + 2*d2Pi  + 2*d3Pi  + d4Pi)
    return new_phi, new_Pi


def evolve_refined_scalar(phi0_c, Pi0_c, phi0_f, Pi0_f,
                          N_c, N_f, offset, dt, n_steps, m2=0.0):
    """Two-grid scalar evolution with simple coupling.

    Args:
        phi0_c, Pi0_c: (N_c, N_c, N_c) coarse initial data
        phi0_f, Pi0_f: (N_f, N_f, N_f) fine initial data
        offset: (ox, oy, oz) coarse index of fine cell (0,0,0)
        dt: timestep (same for both levels; CFL requires dt < h_fine/2)

    Returns:
        dict with final coarse/fine phi and Pi arrays.
    """
    L0 = 1.0  # not used; kept for interface symmetry
    h_c = 1.0  # not used; FD uses unit spacing implicitly in test
    # Use unit spacing for the test — caller provides h via offset scale.
    # In practice grid spacing matters; caller must set h explicitly.
    raise NotImplementedError("use evolve_refined_scalar_with_h instead")


def evolve_refined_scalar_with_h(phi0_c, Pi0_c, phi0_f, Pi0_f,
                                  h_c, h_f, offset, dt, n_steps, m2=0.0):
    """Explicit-spacing version."""
    phi_c, Pi_c = phi0_c.copy(), Pi0_c.copy()
    phi_f, Pi_f = phi0_f.copy(), Pi0_f.copy()
    N_c = phi_c.shape[0]
    N_f = phi_f.shape[0]

    for step in range(n_steps):
        # 1. Evolve coarse
        phi_c, Pi_c = _rk4_step(phi_c, Pi_c, h_c, dt, m2)

        # 2. Prolongate coarse -> fine boundary
        phi_bc_f = prolongate(phi_c, N_f, offset)
        Pi_bc_f = prolongate(Pi_c, N_f, offset)
        # Set fine boundary values (first/last layer)
        for edge in [0, -1]:
            phi_f[edge, :, :] = phi_bc_f[edge, :, :]
            phi_f[:, edge, :] = phi_bc_f[:, edge, :]
            phi_f[:, :, edge] = phi_bc_f[:, :, edge]
            Pi_f[edge, :, :] = Pi_bc_f[edge, :, :]
            Pi_f[:, edge, :] = Pi_bc_f[:, edge, :]
            Pi_f[:, :, edge] = Pi_bc_f[:, :, edge]

        # 3. Evolve fine
        phi_f, Pi_f = _rk4_step(phi_f, Pi_f, h_f, dt, m2)

        # 4. Restrict fine -> coarse (interior)
        phi_c_updated = restrict(phi_f, N_c, offset)
        Pi_c_updated = restrict(Pi_f, N_c, offset)
        # Only overwrite coarse in the region fully covered by fine
        ox, oy, oz = offset
        N_f_half = N_f // 2
        a, b = ox, ox + N_f_half - 1
        # For 2:1 refinement, coarse cell i is fully covered by fine
        # for i in [ox, ox + N_f/2 - 1]
        sl = (slice(a, b+1), slice(a, b+1), slice(a, b+1))
        phi_c[sl] = phi_c_updated[sl]
        Pi_c[sl] = Pi_c_updated[sl]

    return {"phi_coarse": phi_c, "Pi_coarse": Pi_c,
            "phi_fine": phi_f, "Pi_fine": Pi_f}


__all__ = ["evolve_refined_scalar_with_h"]
