"""Binary orbit with complex scalar (v10.0-alpha2).

Complex scalar Psi = Psi_R + i Psi_I. Phase difference between
the two peaks produces anisotropic stress in T_ij, which is
necessary for nonzero TT metric perturbation.

Physics:
  T_ij = 2 Re(d_i Psi* d_j Psi) + gamma_ij L_eff
  With d_i Psi_R perpendicular d_j Psi_I at different points,
  the TT projection is nonzero.
"""
from __future__ import annotations
import numpy as np
from .binary_orbit_ic import circular_orbit_velocity


def complex_binary_ic(N, L, A_R, A_I, sigma, separation, omega, t=0.0):
    """Two complex-scalar Gaussian peaks in circular orbit.

    Body 1: Psi_R peak, Body 2: Psi_I peak. Phase separation
    produces anisotropic stress.
    """
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    r = separation / 2.0
    th = omega * t
    c, s = np.cos(th), np.sin(th)
    x1, y1 = r * c, r * s
    x2, y2 = -r * c, -r * s

    R2_1 = (X - x1)**2 + (Y - y1)**2 + Z**2
    R2_2 = (X - x2)**2 + (Y - y2)**2 + Z**2

    # Body 1 is real-part peak, Body 2 is imaginary-part peak
    PR = A_R * np.exp(-R2_1/(2*sigma**2))
    PI = A_I * np.exp(-R2_2/(2*sigma**2))

    # Momentum Pi_R, Pi_I: tangential motion
    if omega > 0:
        vx1, vy1 = -omega * r * s, omega * r * c
        vx2, vy2 = omega * r * s, -omega * r * c
        dPR_dx = -(X - x1) / sigma**2 * np.exp(-R2_1/(2*sigma**2))
        dPR_dy = -(Y - y1) / sigma**2 * np.exp(-R2_1/(2*sigma**2))
        dPI_dx = -(X - x2) / sigma**2 * np.exp(-R2_2/(2*sigma**2))
        dPI_dy = -(Y - y2) / sigma**2 * np.exp(-R2_2/(2*sigma**2))
        PiR = A_R * (vx1 * dPR_dx + vy1 * dPR_dy)
        PiI = A_I * (vx2 * dPI_dx + vy2 * dPI_dy)
    else:
        PiR = np.zeros_like(PR); PiI = np.zeros_like(PI)

    return {
        "Psi_R": PR, "Psi_I": PI,
        "Pi_R": PiR, "Pi_I": PiI,
        "body1_pos": (x1, y1), "body2_pos": (x2, y2),
        "omega": omega, "separation": separation,
    }


__all__ = ["complex_binary_ic"]
