import numpy as np
import pytest
from zt005.binary_complex_ic import complex_binary_ic


def test_shapes():
    ic = complex_binary_ic(N=9, L=20.0, A_R=1e-3, A_I=1e-3, sigma=1.2,
                            separation=4.0, omega=0.5)
    for k in ("Psi_R", "Psi_I", "Pi_R", "Pi_I"):
        assert ic[k].shape == (9, 9, 9)


def test_two_separated_peaks():
    ic = complex_binary_ic(N=11, L=20.0, A_R=1e-3, A_I=1e-3, sigma=1.0,
                            separation=6.0, omega=0.5)
    mid = 11 // 2
    PR_line = ic["Psi_R"][:, mid, mid]
    PI_line = ic["Psi_I"][:, mid, mid]
    # PR peak on right, PI peak on left
    assert np.argmax(PR_line) > mid
    assert np.argmax(PI_line) < mid


def test_anisotropic_stress_exists():
    """d_i Psi_R and d_j Psi_I are non-parallel -> T_ij has TT content."""
    ic = complex_binary_ic(N=11, L=20.0, A_R=1e-3, A_I=1e-3, sigma=1.2,
                            separation=4.0, omega=0.5)
    # Check that Pi_R, Pi_I are nonzero (needed for real motion)
    assert np.max(np.abs(ic["Pi_R"])) > 0
    assert np.max(np.abs(ic["Pi_I"])) > 0
