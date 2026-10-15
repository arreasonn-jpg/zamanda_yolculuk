"""Brill-Lindquist binary black hole initial data (v11.0-alpha2).

Two punctures at positions +/- d/2 along x-axis, each with mass m.
Time-symmetric (K_ij = 0). Conformal factor:
    psi = 1 + m1/(2 r1) + m2/(2 r2)
This is the Brill-Lindquist (1963) solution of the Hamiltonian
constraint for two non-spinning punctures.

Physical picture:
- Newtonian orbit: f_orbit = (1/2pi) sqrt(G M_total / d^3)
- GW frequency: f_GW = 2 f_orbit
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState

def brill_lindquist_state(N, L, m1=0.5, m2=0.5, d=4.0, eps=None):
    """Brill-Lindquist binary BH IC in Cartesian grid.

    m1, m2: puncture masses
    d: separation along x-axis
    eps: regularization (defaults to dx)
    """
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    dx = L / (N - 1)
    if eps is None:
        eps = dx
    # Puncture positions
    x1 = +d/2.0
    x2 = -d/2.0
    r1 = np.sqrt((X - x1)**2 + Y**2 + Z**2 + eps**2)
    r2 = np.sqrt((X - x2)**2 + Y**2 + Z**2 + eps**2)
    psi = 1.0 + m1/(2.0*r1) + m2/(2.0*r2)
    phi = np.log(psi)
    # Lapse: alpha = (1 - m1/(2r1) - m2/(2r2)) / psi  (approx)
    num = 1.0 - m1/(2.0*r1) - m2/(2.0*r2)
    alpha = np.where(num > 0.01, num / psi, 0.01)
    gtilde = np.broadcast_to(np.eye(3)[None,None,None,:,:],
                             (N, N, N, 3, 3)).copy()
    return BSSNState(
        phi=phi, gtilde=gtilde,
        K=np.zeros((N, N, N)),
        Atilde=np.zeros((N, N, N, 3, 3)),
        Gamma=np.zeros((N, N, N, 3)),
        alpha=alpha,
        beta=np.zeros((N, N, N, 3)),
    )

def orbital_frequency(m1, m2, d):
    """Newtonian orbital frequency (geometrized units G=c=1)."""
    M = m1 + m2
    return (1.0/(2.0*np.pi)) * np.sqrt(M / d**3)

def gw_frequency(m1, m2, d):
    """GW frequency = 2 * orbital frequency (quadrupole)."""
    return 2.0 * orbital_frequency(m1, m2, d)

def adm_mass(m1, m2, d):
    """Brill-Lindquist ADM mass (Newtonian approx)."""
    # Exact: M_ADM = m1 + m2 + interaction correction
    # Leading order: M_ADM = m1 + m2 - m1*m2/d + ...
    return m1 + m2 - m1*m2/d

__all__ = ["brill_lindquist_state", "orbital_frequency",
           "gw_frequency", "adm_mass"]
