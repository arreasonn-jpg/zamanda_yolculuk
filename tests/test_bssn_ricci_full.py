import pytest
pytestmark = pytest.mark.slow
import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_ricci_full import (
    ricci_scalar_physical, hamiltonian_constraint_full,
)
from zt005.bssn_benchmark import schwarzschild_isotropic


def test_minkowski_R_zero():
    s = minkowski_state(7)
    R = ricci_scalar_physical(s.phi, s.gtilde, h=0.1)
    assert np.allclose(R, 0.0, atol=1e-10)


def test_phi_bump_gives_R():
    """Gaussian phi -> R must be nonzero and finite."""
    N = 11
    s = minkowski_state(N)
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.phi = 1e-3 * np.exp(-(X**2+Y**2+Z**2)/0.3)
    R = ricci_scalar_physical(s.phi, s.gtilde, h=2.0/(N-1))
    assert np.all(np.isfinite(R))
    assert np.max(np.abs(R)) > 1e-6, f"R max = {np.max(np.abs(R))}"


def test_schwarzschild_isotropic_R_small():
    """Schwarzschild in isotropic coords is vacuum: R ~ 0."""
    N = 51
    L = 40.0
    h = L / (N - 1)
    r = np.linspace(2.0, 40.0, N)
    M = 1.0
    phi_1d = np.log(1.0 + M/(2.0*r))
    phi = np.broadcast_to(phi_1d[None, None, :], (N, N, N)).copy()
    gt = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    R = ricci_scalar_physical(phi, gt, h)
    interior = R[5:-5, 5:-5, 5:-5]
    assert np.max(np.abs(interior)) < 5e-2, f"max|R|={np.max(np.abs(interior))}"


def test_hamiltonian_full_shapes():
    s = minkowski_state(7)
    H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h=0.1)
    assert H.shape == (7,7,7)
