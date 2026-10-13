import numpy as np
import pytest
from zt005.gw_extraction import (
    metric_perturbation_from_state, extract_at_observer,
    extract_multiple_observers, psi4_from_history,
    strain_amplitude_from_psi4_history,
)
from zt005.bssn_state import minkowski_state


def test_metric_perturbation_zero_for_minkowski():
    s = minkowski_state(5)
    h = metric_perturbation_from_state(s)
    assert np.allclose(h, 0.0, atol=1e-12)


def test_extract_zero_for_minkowski():
    s = minkowski_state(5)
    k = np.array([0.0, 0.0, 1.0])
    e1 = np.array([1.0, 0.0, 0.0])
    e2 = np.array([0.0, 1.0, 0.0])
    hp, hc = extract_at_observer(s, (2, 2, 2), k, e1, e2)
    assert hp == pytest.approx(0.0, abs=1e-12)
    assert hc == pytest.approx(0.0, abs=1e-12)


def test_psi4_from_history_sinusoidal():
    """h_+ = A sin(omega t) -> psi_4 = -1/2 * h_ddot."""
    A = 1e-20
    omega = 1e3
    dt = 1e-5
    t = np.arange(0, 1e-3, dt)
    hp = A * np.sin(omega * t)
    hc = np.zeros_like(t)
    psi4 = psi4_from_history(hp, hc, dt)
    assert len(psi4) > 0
    # |psi_4| ~ 1/2 omega^2 A
    expected = 0.5 * omega**2 * A
    assert np.max(np.abs(psi4)) == pytest.approx(expected, rel=0.1)


def test_multiple_observers():
    s = minkowski_state(5)
    k = np.array([0.0, 0.0, 1.0])
    e1 = np.array([1.0, 0.0, 0.0])
    e2 = np.array([0.0, 1.0, 0.0])
    obs = [(1,1,1), (2,2,2), (3,3,3)]
    out = extract_multiple_observers(s, obs, k, e1, e2)
    assert len(out) == 3
    for entry in out:
        assert "h_plus" in entry
        assert "h_cross" in entry


def test_strain_recovery_from_psi4():
    """Sinusoidal psi_4 -> strain ~ A."""
    A = 1e-20
    omega = 1e3
    dt = 1e-5
    t = np.arange(0, 1e-2, dt)
    psi4 = 0.5 * omega**2 * A * np.exp(-1j * omega * t)
    h_est = strain_amplitude_from_psi4_history(psi4, dt)
    assert len(h_est) > 0
    assert np.median(h_est[100:-100]) == pytest.approx(A, rel=0.2)
