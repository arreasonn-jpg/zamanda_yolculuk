import numpy as np
from zt005.matter.matter_brans_dicke import brans_dicke_sources, gr_limit_check
from zt005.bssn_state import minkowski_state


def test_vacuum_zero():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    src = brans_dicke_sources(Z + 1.0, Z, s.gtilde, s.phi, h=0.1)
    assert np.allclose(src["rho"], 0.0, atol=1e-12)
    assert np.allclose(src["S_i"], 0.0, atol=1e-12)


def test_constant_phi_gives_zero():
    """phi = const, Pi = 0 -> rho = 0 (only G_eff modified)."""
    s = minkowski_state(5)
    phi_BD = np.ones((5,5,5)) * 2.0
    Pi_BD = np.zeros((5,5,5))
    src = brans_dicke_sources(phi_BD, Pi_BD, s.gtilde, s.phi, h=0.1)
    assert np.allclose(src["rho"], 0.0, atol=1e-12)


def test_pi_only_gives_positive_rho():
    s = minkowski_state(5)
    phi_BD = np.ones((5,5,5))
    Pi_BD = np.ones((5,5,5)) * 0.5
    src = brans_dicke_sources(phi_BD, Pi_BD, s.gtilde, s.phi,
                                h=0.1, omega=100.0)
    # rho = (omega/phi) Pi^2 = 100 * 0.25 = 25
    assert np.allclose(src["rho"], 25.0)
    assert np.allclose(src["G_eff_factor"], 1.0)


def test_phi_gradient_gives_rho():
    N = 5
    s = minkowski_state(N)
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    phi_BD = 1.0 + X * 0.1
    Pi_BD = np.zeros((N,N,N))
    src = brans_dicke_sources(phi_BD, Pi_BD, s.gtilde, s.phi,
                                h=0.5, omega=100.0)
    # interior has grad2 = 0.01, rho = 100/1 * 0.01 = 1
    interior = src["rho"][2:-2, 2:-2, 2:-2]
    assert np.allclose(interior, 1.0, atol=0.5)


def test_gr_limit():
    assert gr_limit_check(50000.0)
    assert not gr_limit_check(1000.0)


def test_shapes():
    s = minkowski_state(7)
    Z = np.zeros((7,7,7))
    src = brans_dicke_sources(Z + 1.0, Z, s.gtilde, s.phi, h=0.1)
    assert src["rho"].shape == (7,7,7)
    assert src["S_i"].shape == (7,7,7,3)
    assert src["S_ij"].shape == (7,7,7,3,3)
