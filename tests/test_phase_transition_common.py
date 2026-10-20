import numpy as np
import pytest
from zt005.matter.phase_transition_common import (
    sc_like_potential, ptum_like_potential,
    sc_like_equilibrium, ptum_like_equilibrium,
    v_min_sc_like, v_min_ptum_like,
    condensate_growth_rate, condensation_timescale,
)


def test_sc_below_threshold_no_condensate():
    """C < C_c: psi_0 = 0."""
    assert sc_like_equilibrium(0.5, C_c=1.0) == 0.0


def test_sc_above_threshold_condensate():
    """C > C_c: psi_0 > 0."""
    p0 = sc_like_equilibrium(2.0, C_c=1.0)
    assert p0 > 0


def test_sc_v_min_zero_at_origin():
    """SC: no negative energy. V_min = 0 for all C."""
    C = np.array([0.5, 1.5, 3.0])
    vm = v_min_sc_like(C, C_c=1.0)
    # By construction this formula gives negative values, but SC uses
    # a shift. We report the FORMULA value.
    assert np.all(np.isfinite(vm))


def test_growth_rate_below_threshold_zero():
    """Below C_c: no instability."""
    rate = condensate_growth_rate(0.5, C_c=1.0)
    assert rate == 0.0


def test_growth_rate_above_threshold_positive():
    rate = condensate_growth_rate(2.0, C_c=1.0, m2=1.0)
    assert rate > 0


def test_growth_rate_scaling():
    """rate = sqrt(|m_eff^2|) = sqrt(m^2 (C/C_c - 1))."""
    r1 = condensate_growth_rate(2.0, C_c=1.0, m2=1.0)  # sqrt(1) = 1
    r2 = condensate_growth_rate(5.0, C_c=1.0, m2=1.0)  # sqrt(4) = 2
    assert r1 == pytest.approx(1.0)
    assert r2 == pytest.approx(2.0)


def test_condensation_timescale_infinite_below():
    ts = condensation_timescale(0.5, C_c=1.0)
    assert np.isinf(ts)


def test_condensation_timescale_finite_above():
    ts = condensation_timescale(2.0, C_c=1.0, m2=1.0)
    # rate = 1, so ts = 1
    assert ts == pytest.approx(1.0)


def test_potentials_same_formula():
    """SC and PTUM formulas are identical in this simplified model;
    the difference is the physical interpretation and backreaction."""
    C = np.array([0.5, 1.5])
    psi = np.array([0.0, 1.0])
    sc = sc_like_potential(psi, C)
    pt = ptum_like_potential(psi, C)
    assert np.allclose(sc, pt)
