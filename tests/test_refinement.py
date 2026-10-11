import numpy as np
import pytest
from zt005.refinement import prolongate, restrict


def test_prolongate_constant():
    N0 = 11
    f = np.ones((N0, N0, N0)) * 5.0
    g = prolongate(f, 5, (3, 3, 3))
    assert np.allclose(g, 5.0, atol=1e-12)


def test_prolongate_linear_exact():
    """Linear f = 2X + 1: g[i,j,k] should = 2*(3+i/2)+1 = 7+i."""
    N0 = 11
    x = np.arange(N0, dtype=float)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    f = 2.0 * X + 1.0
    g = prolongate(f, 5, (3, 3, 3))
    for i in range(5):
        expected = 7.0 + i
        assert np.allclose(g[i, :, :], expected, atol=1e-12)


def test_restrict_constant():
    N_f = 5
    f = np.ones((N_f, N_f, N_f)) * 3.0
    g = restrict(f, 11, (3, 3, 3))
    assert g[4, 4, 4] == pytest.approx(3.0, abs=1e-12)


def test_prolongate_shapes():
    f = np.ones((21, 21, 21))
    g = prolongate(f, 5, (8, 8, 8))
    assert g.shape == (5, 5, 5)


def test_restrict_shapes():
    f = np.ones((5, 5, 5))
    g = restrict(f, 21, (8, 8, 8))
    assert g.shape == (21, 21, 21)


def test_prolongate_smooth_roundtrip():
    """Smooth f: prolongate + restrict ~= identity in fine region."""
    N0 = 21
    x = np.linspace(0, 1, N0)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    f = np.sin(3*X) * np.cos(2*Y) + 0.5*Z
    offset = (8, 8, 8)
    g = prolongate(f, 5, offset)
    f_back = restrict(g, N0, offset)
    # cells [8..10] fully covered by fine
    a = offset[0]
    b = offset[0] + 3
    diff = np.abs(f_back[a:b, a:b, a:b] - f[a:b, a:b, a:b])
    assert diff.max() < 0.5, f"max diff = {diff.max()}"
