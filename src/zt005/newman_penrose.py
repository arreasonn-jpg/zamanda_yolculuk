"""Newman-Penrose scalars for GW extraction (v9.0-alpha1).

In linearized GR with TT-gauge metric perturbation:
    h_ij = h_+ e_+_ij + h_x e_x_ij

The NP scalar psi_4 is:
    psi_4 = -1/2 (h_+ddot - i h_xddot)

For a monochromatic wave with h ~ A exp(-i omega t):
    |psi_4| = 1/2 omega^2 |h|

We implement:
- TT-gauge extraction of h_+, h_x from a BSSN metric perturbation
- psi_4 from the second time derivative
- GW energy flux (Peters-Mathews formula)
"""
from __future__ import annotations
import numpy as np

G = 6.67430e-11
C = 2.99792458e8


def tt_polarizations(h_tensor, k_hat, e1, e2):
    """Extract h_+ and h_x from a metric perturbation tensor.

    h_tensor: 3x3 symmetric tensor (spatial part of h_munu)
    k_hat: propagation direction (unit vector)
    e1, e2: transverse orthonormal vectors (e1 perpendicular e2, both perp k_hat)

    h_+ = e1^T h e1 - e2^T h e2
    h_x = e1^T h e2 + e2^T h e1
    """
    h = np.asarray(h_tensor, float)
    e1 = np.asarray(e1, float)
    e2 = np.asarray(e2, float)
    h_plus = float(e1 @ h @ e1 - e2 @ h @ e2)
    h_cross = float(e1 @ h @ e2 + e2 @ h @ e1)
    return h_plus, h_cross


def psi4_from_hddot(h_plus_ddot, h_cross_ddot):
    """psi_4 = -1/2 (h_+ddot - i h_xddot)."""
    return -0.5 * (h_plus_ddot - 1j * h_cross_ddot)


def psi4_amplitude(h_amp, omega):
    """|psi_4| = 1/2 omega^2 |h| for a monochromatic wave."""
    return 0.5 * omega**2 * abs(h_amp)


def gw_energy_flux(h_plus_dot, h_cross_dot):
    """GW energy flux (Peters-Mathews):
    F = (c^3 / (16 pi G)) < (h_+dot)^2 + (h_xdot)^2 >
    """
    return (C**3 / (16.0 * np.pi * G)) * (h_plus_dot**2 + h_cross_dot**2)


def strain_from_psi4(psi4, omega):
    """Invert |psi_4| = 1/2 omega^2 |h| to get |h|."""
    return 2.0 * abs(psi4) / max(omega**2, 1e-30)


__all__ = [
    "tt_polarizations", "psi4_from_hddot", "psi4_amplitude",
    "gw_energy_flux", "strain_from_psi4",
]
