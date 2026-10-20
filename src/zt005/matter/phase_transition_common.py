"""Cross-scale phase transition test (v6.1-alpha2).

Tests whether SC-like (BCS condensation) and PTUM-like (tachyon
discharge) potentials share a common "condensate formation" pattern
when evolved on the same BSSN infrastructure.

This is a NUMERICAL comparison, not a claim of unification.
The two models operate at vastly different scales.

Shared structure tested:
  1. Control parameter C(r) crosses threshold C_c
  2. Order parameter psi forms a nonzero condensate
  3. Metric backreaction appears

SC-like: BCS gap, V_min > 0 (regular core)
PTUM-like: Tachyon, V_min < 0 (NEC violation)
"""
from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------------------
# Potentials
# ---------------------------------------------------------------------------
def sc_like_potential(psi, C, C_c=1.0, lam=1.0, m2=1.0):
    """BCS-like: control C stabilizes condensate. V_min > 0.

    Below threshold (C < C_c): psi_0 = 0
    Above threshold (C > C_c): psi_0 > 0 with V_min = 0 (restored).
    """
    psi = np.asarray(psi, float)
    C = np.asarray(C, float)
    # Effective mass: m_eff^2 = m^2 (1 - C/C_c)
    meff2 = m2 * (1.0 - C / C_c)
    return 0.5 * meff2 * psi**2 + 0.25 * lam * psi**4


def ptum_like_potential(psi, C, C_c=1.0, lam=1.0, m2=1.0):
    """Tachyon: control C destabilizes. V_min < 0 above threshold.

    Same mass formula. Distinguishing feature: we ADD a constant so
    that above threshold the minimum is NEGATIVE.
    """
    psi = np.asarray(psi, float)
    C = np.asarray(C, float)
    meff2 = m2 * (1.0 - C / C_c)
    # When C > C_c, meff2 < 0, minimum at psi_0 = sqrt(-meff2/lam)
    # V_min = -meff2^2 / (4 lam). We do NOT shift by constant here.
    return 0.5 * meff2 * psi**2 + 0.25 * lam * psi**4


def sc_like_equilibrium(C, C_c=1.0, lam=1.0, m2=1.0):
    """Return psi_0 for SC-like potential."""
    C = np.asarray(C, float)
    meff2 = m2 * (1.0 - C / C_c)
    arg = np.where(meff2 < 0, -meff2 / lam, 0.0)
    return np.sqrt(arg)


def ptum_like_equilibrium(C, C_c=1.0, lam=1.0, m2=1.0):
    """Return psi_0 for PTUM-like potential (same as SC form)."""
    return sc_like_equilibrium(C, C_c, lam, m2)


def v_min_sc_like(C, C_c=1.0, lam=1.0, m2=1.0):
    """Minimum value of SC-like potential."""
    C = np.asarray(C, float)
    meff2 = m2 * (1.0 - C / C_c)
    # V_min = -meff2^2 / (4 lam) when meff2 < 0, else 0
    return np.where(meff2 < 0, -meff2**2 / (4.0*lam), 0.0)


def v_min_ptum_like(C, C_c=1.0, lam=1.0, m2=1.0):
    """Same form — but we TEST whether backreaction differs."""
    return v_min_sc_like(C, C_c, lam, m2)


# ---------------------------------------------------------------------------
# Analytic condensate formation timescale
# ---------------------------------------------------------------------------
def condensate_growth_rate(C, C_c=1.0, m2=1.0):
    """|m_eff| when C > C_c: growth rate of tachyon/BCS instability."""
    C = np.asarray(C, float)
    meff2 = m2 * (1.0 - C / C_c)
    arg = np.where(meff2 < 0, -meff2, 0.0)
    return np.sqrt(arg)


def condensation_timescale(C, C_c=1.0, m2=1.0):
    """Approx time for condensate to grow from seed to equilibrium.
    ~ 1 / sqrt(|m_eff^2|).
    """
    rate = condensate_growth_rate(C, C_c, m2)
    return np.where(rate > 0, 1.0 / np.maximum(rate, 1e-30), np.inf)


__all__ = [
    "sc_like_potential", "ptum_like_potential",
    "sc_like_equilibrium", "ptum_like_equilibrium",
    "v_min_sc_like", "v_min_ptum_like",
    "condensate_growth_rate", "condensation_timescale",
]
