"""ADM mass extraction for asymptotically flat BSSN slices.

  M_ADM = (1/(16 pi)) oint (d_i g_jj - d_j g_ij) n^i dS

We use the simpler volume-integrated form (Baumgarte-Shapiro 2010, eq 2.135):

  M_ADM = (1/(16 pi)) int [ d_i d_j g_ij - d^2 g_ii ] dV

In BSSN variables (approx for slowly-varying psi):
  M_ADM ~ -(1/(2 pi)) int d^2 psi dV
where psi = e^phi.
"""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central


def adm_mass_approx(phi, h):
    """Approximate ADM mass from conformal factor phi.

    For Schwarzschild isotropic: psi = 1 + M/(2r), laplacian(psi) = 0
    except at r=0 (delta). Volume integral of -laplacian(psi)/(2 pi) = M.
    Numerically we compute laplacian(psi) and integrate over interior.
    """
    psi = np.exp(phi)
    d = d_central(psi, h)
    lap = np.zeros_like(psi)
    for ax in range(3):
        dd = d_central(d[ax], h)
        lap += dd[ax]
    Nx, Ny, Nz = psi.shape
    M = -(1.0 / (2.0 * np.pi)) * np.sum(lap) * h**3
    return float(M)


__all__ = ["adm_mass_approx"]
