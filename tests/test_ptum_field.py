import numpy as np
import pytest
from zt005.ptum_field import (
    V_fuel, dV_fuel, m_eff_squared, Phi_critical,
    V_proxy, psi_equilibrium, fuel_profile, ptdm_equilibrium_check,
)


def test_V_fuel_zero_at_origin():
    assert V_fuel(0.0) == 0.0


def test_V_fuel_plateau():
    """At large Phi, V -> V0."""
    V_large = V_fuel(10.0, V0=1.0, lam=0.05)
    assert V_large == pytest.approx(1.0, abs=1e-2)


def test_dV_fuel_positive_for_positive_Phi():
    assert dV_fuel(1.0) > 0
    assert dV_fuel(-1.0) < 0


def test_m_eff_squared_sign_change():
    """Sign flips at Phi_c."""
    Phi_c = Phi_critical(alpha=1.0, Eth=1.0)
    assert m_eff_squared(0.5, alpha=1.0, Eth=1.0) > 0  # dormant
    assert m_eff_squared(2.0, alpha=1.0, Eth=1.0) < 0  # activated


def test_psi_equilibrium_zero_when_dormant():
    """Below threshold, psi_0 = 0."""
    assert psi_equilibrium(0.5, alpha=1.0, Eth=1.0) == 0.0


def test_psi_equilibrium_nonzero_when_activated():
    """Above threshold, psi_0 > 0."""
    psi0 = psi_equilibrium(2.0, m2=1.0, alpha=1.0, Eth=1.0, lam_psi=1.0)
    assert psi0 > 0


def test_fuel_profile_gaussian_shape():
    phi = fuel_profile(N=11, L=20.0, Phi0=0.5, sigma=2.0)
    # center value = Phi0
    assert phi[5,5,5] == pytest.approx(0.5, rel=1e-6)
    # edges < center
    assert phi[0,0,0] < phi[5,5,5]


def test_ptdm_check_matches_threshold():
    Phi = np.array([0.5, 1.0, 2.0])
    act = ptdm_equilibrium_check(Phi, alpha=1.0, Eth=1.0)
    # Phi=1.0 exactly at threshold -> m_eff=0, not < 0
    assert act.tolist() == [False, False, True]
