import numpy as np
import pytest
from zt005.bssn_benchmark import schwarzschild_isotropic, ricci_scalar_radial
from zt005.bssn import ADM, adm_to_bssn


def test_schwarzschild_gtilde_flat():
    s = schwarzschild_isotropic(M=1.0, r=10.0)
    assert np.allclose(s["gtilde"], np.eye(3))
    assert s["K"] == 0.0
    assert np.allclose(s["Atilde"], 0.0)


def test_schwarzschild_phi_positive():
    s = schwarzschild_isotropic(M=1.0, r=10.0)
    assert s["phi"] > 0.0
    # psi = 1 + M/(2r) = 1.05 -> phi = ln(1.05)
    assert s["phi"] == pytest.approx(np.log(1.05), rel=1e-12)


def test_schwarzschild_alpha_limits():
    # far field: alpha -> 1
    s = schwarzschild_isotropic(M=1.0, r=1e6)
    assert s["alpha"] == pytest.approx(1.0, abs=1e-5)
    # near horizon: alpha -> 0
    s2 = schwarzschild_isotropic(M=1.0, r=0.5001)
    assert s2["alpha"] < 1e-3


def test_minkowski_ricci_zero():
    r = np.linspace(1.0, 5.0, 201)
    h = r[1] - r[0]
    phi = np.zeros_like(r)
    R = ricci_scalar_radial(phi, r, h)
    assert np.allclose(R[5:-5], 0.0, atol=1e-12)


def test_schwarzschild_ricci_small():
    """Vacuum Schwarzschild: R ~ 0 to FD accuracy."""
    r = np.linspace(2.0, 20.0, 2001)
    h = r[1] - r[0]
    M = 1.0
    phi = np.log(1.0 + M / (2.0 * r))
    R = ricci_scalar_radial(phi, r, h)
    # Check interior points away from boundaries
    interior = R[50:-50]
    assert np.max(np.abs(interior)) < 1e-3
