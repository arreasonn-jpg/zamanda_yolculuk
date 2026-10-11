import numpy as np
import pytest
from zt005.kerr_sc_ic import kerr_sc_metric, kerr_sc_state


def test_a0_center_regular():
    """M_eff(0) = 0 -> H(0) = 0 -> alpha(0) = 1."""
    met = kerr_sc_metric(N=9, L=20.0, M=1.0, a=0.0, rc=2.0)
    # center index (4,4,4)
    assert met["alpha"][4,4,4] == pytest.approx(1.0, abs=1e-6)
    assert met["H"][4,4,4] < 1e-6


def test_far_field_recovers_kerr():
    """At r >> rc, M_eff -> M and metric matches standard KS."""
    met = kerr_sc_metric(N=9, L=20.0, M=1.0, a=0.5, rc=0.5)
    # corner of grid: r ~ sqrt(3)*10
    # M_eff should be ~1 there
    assert met["M_eff"][0,0,0] == pytest.approx(1.0, rel=1e-3)


def test_a0_beta_is_radial():
    """a=0: KS shift is radial (nonzero), no azimuthal component.

    beta_i should be parallel to (x, y, z) direction at each point.
    """
    met = kerr_sc_metric(N=9, L=20.0, M=1.0, a=0.0, rc=2.0)
    # beta_i nonzero (KS radial shift) but not azimuthal
    assert np.max(np.abs(met["beta_i"])) > 1e-6
    # Cross product beta x r should be zero for pure radial shift
    x = np.linspace(-10, 10, 9)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.stack([X, Y, Z], axis=-1)
    cross = np.cross(met["beta_i"], r)
    interior = cross[2:-2, 2:-2, 2:-2]
    assert np.allclose(interior, 0.0, atol=1e-10)


def test_a_nonzero_beta():
    met = kerr_sc_metric(N=9, L=20.0, M=1.0, a=0.9, rc=2.0)
    assert np.max(np.abs(met["beta_i"])) > 1e-6


def test_alpha_bounded():
    met = kerr_sc_metric(N=9, L=20.0, M=1.0, a=0.9, rc=2.0)
    assert np.all(met["alpha"] > 0)
    assert np.all(met["alpha"] <= 1.0 + 1e-10)


def test_state_returns_bssn():
    s = kerr_sc_state(N=9, L=20.0, M=1.0, a=0.5, rc=2.0)
    assert s.phi.shape == (9,9,9)
    assert s.beta.shape == (9,9,9,3)
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(s.K))


def test_a_too_large_raises():
    with pytest.raises(ValueError):
        kerr_sc_state(N=9, L=20.0, M=1.0, a=1.5, rc=2.0)


def test_det_gtilde_unity():
    s = kerr_sc_state(N=9, L=20.0, M=1.0, a=0.5, rc=2.0)
    dets = np.linalg.det(s.gtilde)
    assert np.allclose(dets, 1.0, atol=1e-10)


def test_rc_dependence():
    """Smaller rc: sharper core transition."""
    met_small = kerr_sc_metric(N=9, L=20.0, M=1.0, a=0.5, rc=0.5)
    met_large = kerr_sc_metric(N=9, L=20.0, M=1.0, a=0.5, rc=3.0)
    # Small rc: at corner M_eff already close to M
    # Large rc: at corner M_eff < M
    assert met_small["M_eff"][0,0,0] > met_large["M_eff"][0,0,0]
