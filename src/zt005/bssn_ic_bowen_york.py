"""Bowen-York binary BH initial data (v12.0-alpha1).

Two punctures with momentum P (tangential) -> orbital motion.
For equal-mass circular orbit (Newtonian approx):
    P = m * v_orbit,  v_orbit = sqrt(M/(4d))  (M = m1+m2)
The Bowen-York conformal extrinsic curvature:
    A_ij = (3/(2 r^2)) [ P_i n_j + P_j n_i - (delta_ij - n_i n_j) P_k n_k ]
    (simplified, one puncture, flat space)

For two punctures with opposite momenta, superpose.
Note: this is the flat-space approximation. Full Bowen-York
solution requires solving the momentum constraint; we use the
leading-order superposition as IC.
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState

def bowen_york_binary(N, L, m1=0.5, m2=0.5, d=6.0, eps=None):
    """Two equal-mass punctures with tangential momenta.
    Direction: counter-clockwise in xy plane.
    """
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    dx = L / (N - 1)
    if eps is None:
        eps = dx
    x1 = +d/2.0; x2 = -d/2.0
    r1 = np.sqrt((X - x1)**2 + Y**2 + Z**2 + eps**2)
    r2 = np.sqrt((X - x2)**2 + Y**2 + Z**2 + eps**2)
    # Conformal factor (same as BL)
    psi = 1.0 + m1/(2.0*r1) + m2/(2.0*r2)
    phi = np.log(psi)
    # Orbital velocity (Newtonian, geometrized units)
    M_total = m1 + m2
    v_orb = np.sqrt(M_total / (4.0 * d))
    # Momenta: tangential (y-direction)
    P1 = np.array([0.0, +m1*v_orb, 0.0])
    P2 = np.array([0.0, -m2*v_orb, 0.0])
    # Lapse
    num = 1.0 - m1/(2.0*r1) - m2/(2.0*r2)
    alpha = np.where(num > 0.01, num / psi, 0.01)
    # gtilde = I (conformal flat, leading order)
    gtilde = np.broadcast_to(np.eye(3)[None,None,None,:,:],
                             (N, N, N, 3, 3)).copy()
    # Atilde from Bowen-York formula
    def by_tensor(P, Xi, Yi, Zi, xi, ri_soft):
        # Proper unit vector: use hard radius (no eps) for direction
        ri_hard = np.sqrt((Xi - xi)**2 + Yi**2 + Zi**2 + 1e-20)
        nx = (Xi - xi) / ri_hard
        ny = Yi / ri_hard
        nz = Zi / ri_hard
        Pn = P[0]*nx + P[1]*ny + P[2]*nz
        A = np.zeros((N, N, N, 3, 3))
        n = [nx, ny, nz]
        for i in range(3):
            for j in range(3):
                term1 = P[i]*n[j] + P[j]*n[i]
                delta_ij = 1.0 if i == j else 0.0
                term2 = (delta_ij - n[i]*n[j]) * Pn
                A[..., i, j] = (3.0/(2.0*ri_soft**2)) * (term1 - term2)
        return A
    A1 = by_tensor(P1, X, Y, Z, x1, r1)
    A2 = by_tensor(P2, X, Y, Z, x2, r2)
    # Atilde = psi^{-4} * A_total  (conformal rescaling)
    Atilde = (A1 + A2) * np.exp(-4.0 * phi[..., None, None])
    return BSSNState(
        phi=phi, gtilde=gtilde,
        K=np.zeros((N, N, N)),
        Atilde=Atilde,
        Gamma=np.zeros((N, N, N, 3)),
        alpha=alpha,
        beta=np.zeros((N, N, N, 3)),
    )

__all__ = ["bowen_york_binary"]
