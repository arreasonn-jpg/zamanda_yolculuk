"""PTUM field content for BSSN (v5.0-alpha1).

Reference: Esa (2026), "Phase Transition Universes Model",
eqs (3)-(10), (17).

Minimal field set:
  Phi (fuel)      : plateau potential V_eff(Phi) = V0 (1 - e^{-lambda Phi^2})
  psi (proxy)     : effective mass m_eff^2(Phi) = m^2 (1 - alpha Phi^2 / Eth^2)
                    activated (tachyonic) when Phi > Phi_c = Eth / sqrt(alpha)
  Phi_res (DE)    : W(Phi_res) = W0 e^{-kappa Phi_res^2}(1 - e^{-lambda_res Phi_res^2})

We implement only the fuel + proxy sector for BSSN coupling. The
residual dark-energy field is left as a background (v5.0-alpha2).
"""
from __future__ import annotations
import numpy as np

C = 2.99792458e8


def V_fuel(Phi, V0=1.0, lam=0.05):
    """V_eff(Phi) = V0 (1 - e^{-lambda Phi^2})."""
    Phi = np.asarray(Phi, float)
    return V0 * (1.0 - np.exp(-lam * Phi**2))


def dV_fuel(Phi, V0=1.0, lam=0.05):
    """dV/dPhi = 2 V0 lambda Phi e^{-lambda Phi^2}."""
    Phi = np.asarray(Phi, float)
    return 2.0 * V0 * lam * Phi * np.exp(-lam * Phi**2)


def m_eff_squared(Phi, m2=1.0, alpha=1.0, Eth=1.0):
    """Effective mass squared of proxy field: m^2 (1 - alpha Phi^2/Eth^2)."""
    Phi = np.asarray(Phi, float)
    return m2 * (1.0 - alpha * Phi**2 / Eth**2)


def Phi_critical(alpha=1.0, Eth=1.0):
    """Tachyonic threshold: Phi_c = Eth / sqrt(alpha)."""
    return Eth / np.sqrt(alpha)


def V_proxy(psi, Phi, m2=1.0, alpha=1.0, Eth=1.0, lam_psi=1.0):
    """V_act(psi; Phi) = 1/2 m_eff^2(Phi) psi^2 + 1/4 lam_psi psi^4."""
    psi = np.asarray(psi, float)
    meff2 = m_eff_squared(Phi, m2, alpha, Eth)
    return 0.5 * meff2 * psi**2 + 0.25 * lam_psi * psi**4


def psi_equilibrium(Phi, m2=1.0, alpha=1.0, Eth=1.0, lam_psi=1.0):
    """Post-activation minimum psi_0^2 = -m_eff^2 / lam_psi (0 if dormant)."""
    meff2 = m_eff_squared(Phi, m2, alpha, Eth)
    arg = np.where(meff2 < 0.0, -meff2 / lam_psi, 0.0)
    return np.sqrt(arg)


def fuel_profile(N, L, Phi0=0.5, sigma=2.0):
    """Gaussian fuel profile centered at origin."""
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2
    return Phi0 * np.exp(-R2 / (2.0 * sigma**2))


def ptdm_equilibrium_check(Phi, m2=1.0, alpha=1.0, Eth=1.0, lam_psi=1.0):
    """Return True where tachyonic (activated), False where dormant."""
    return m_eff_squared(Phi, m2, alpha, Eth) < 0.0


__all__ = [
    "V_fuel", "dV_fuel", "m_eff_squared", "Phi_critical",
    "V_proxy", "psi_equilibrium", "fuel_profile", "ptdm_equilibrium_check",
]


def psi_equilibrium_smooth(Phi, m2=1.0, alpha=1.0, Eth=1.0, lam_psi=1.0,
                            width=0.3):
    """Smooth version of psi_equilibrium for use in evolution.

    The sharp discontinuity in psi_equilibrium at Phi = Phi_c creates
    a Laplacian spike that drives numerical instability. We use a
    smoothstep transition of width `width` around Phi_c to regularize.
    """
    Phi = np.asarray(Phi, float)
    meff2 = m2 * (1.0 - alpha * Phi**2 / Eth**2)
    arg = np.where(meff2 < 0.0, -meff2 / lam_psi, 0.0)
    psi_eq = np.sqrt(arg)
    # Smooth transition: 0 for Phi < Phi_c - width, 1 for Phi > Phi_c + width
    # smoothstep = 3 t^2 - 2 t^3
    Phi_c = Eth / np.sqrt(alpha)
    t = (Phi - (Phi_c - width)) / (2.0 * width)
    t = np.clip(t, 0.0, 1.0)
    smooth = 3.0 * t**2 - 2.0 * t**3
    return psi_eq * smooth


__all__ = list(__all__) + ["psi_equilibrium_smooth"]
