"""SC condensate as matter source for BSSN (v4.0-alpha2).

Reference: Esa (2026), Schwarzschild Condensation, eqs (9), (16):
  L_SC = sqrt(-g) [g^mu_nu d_mu Psi* d_nu Psi - m_pair^2 |Psi|^2
                   - (1/2) g_int |Psi|^4]
  Gaussian ground state in de Sitter trap:
    Psi_0(r) = (1/(pi^{3/2} sigma^3))^{1/2} exp(-r^2/(2 sigma^2))
    sigma = sqrt(hbar / (m_pair omega))
    omega = c / ell,   ell = de Sitter length

We use a real scalar approximation (Pi = 0 initially):
  phi_field(r) = A_amp * exp(-r^2 / (2 sigma^2))
with A_amp chosen so that rho_core = (1/2) m_pair^2 phi^2
matches the SC core density.
"""
from __future__ import annotations
import numpy as np

C = 2.99792458e8
HBAR = 1.054571817e-34
G = 6.67430e-11


def sc_condensate_profile(N, L, M=1.0, rho_core=None, m_pair_MeV=10.0,
                          ell=None):
    """Return (phi_field, Pi_field, sigma, omega, A_amp)."""
    # de Sitter length from rho_core
    if rho_core is None:
        # default: use SC core params (from sc_core.py)
        from .sc_core import core_parameters
        p = core_parameters(M, 1e19)
        ell = p["ell"]
        rho_core = 1e19
    elif ell is None:
        ell = np.sqrt(3.0 * C**4 / (8.0 * np.pi * G * rho_core))

    omega = C / ell
    m_pair = m_pair_MeV * 1e6 * 1.602e-19 / C**2  # kg
    sigma = np.sqrt(HBAR / (m_pair * omega))

    # Amplitude: rho_core = (1/2) m_pair^2 phi^2 (leading term)
    A_amp = np.sqrt(2.0 * rho_core) / max(m_pair, 1e-30)

    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2

    # Gaussian: A_amp * exp(-r^2 / (2 sigma^2))
    # For numerical grid, use a smoothed profile in grid units
    sigma_grid = max(sigma / (L / (N - 1)), 0.5)  # at least half a cell
    phi_field = A_amp * np.exp(-R2 / (2.0 * (sigma_grid * L / (N-1))**2))
    Pi_field = np.zeros_like(phi_field)

    return phi_field, Pi_field, sigma, omega, A_amp


def sc_condensate_potential(m_pair_MeV=10.0, g_int=0.1):
    """Return V(phi) = (1/2) m^2 phi^2 + (1/4) g_int phi^4."""
    m_pair = m_pair_MeV * 1e6 * 1.602e-19 / C**2
    m2 = m_pair**2
    def V(p):
        return 0.5 * m2 * p**2 + 0.25 * g_int * p**4
    return V


__all__ = ["sc_condensate_profile", "sc_condensate_potential"]
