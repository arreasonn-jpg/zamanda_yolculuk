import numpy as np
import pytest
from zt005.loop_quantum_cosmology import (
    RHO_PLANCK, RHO_C, friedmann_lqc, friedmann_classical,
    bounce_density, is_bouncing, hubble_ratio,
)


def test_planck_density_order():
    """rho_Planck ~ 5.16e96 kg/m^3."""
    assert 1e96 < RHO_PLANCK < 1e97


def test_lqc_rho_c_order():
    """rho_c ~ 0.41 rho_Planck."""
    assert RHO_C == pytest.approx(0.41 * RHO_PLANCK)


def test_classical_limit_low_rho():
    """rho << rho_c: LQC ~ classical."""
    rho = 1e-10 * RHO_C
    H2_lqc = friedmann_lqc(rho)
    H2_cl = friedmann_classical(rho)
    assert H2_lqc == pytest.approx(H2_cl, rel=1e-9)


def test_bounce_at_rho_c():
    """H^2 = 0 at rho = rho_c."""
    assert friedmann_lqc(RHO_C) == pytest.approx(0.0, abs=1e-30)


def test_bounce_density():
    assert bounce_density() == RHO_C


def test_is_bouncing():
    assert not is_bouncing(0.3 * RHO_C)
    assert is_bouncing(0.6 * RHO_C)
    assert is_bouncing(RHO_C)


def test_hubble_ratio():
    assert hubble_ratio(0.0) == 1.0
    assert hubble_ratio(0.5 * RHO_C) == pytest.approx(0.5)
    assert hubble_ratio(RHO_C) == pytest.approx(0.0)


def test_lqc_negative_for_rho_above_rho_c():
    """Formally H^2 < 0 above rho_c: quantum regime."""
    assert friedmann_lqc(1.5 * RHO_C) < 0
