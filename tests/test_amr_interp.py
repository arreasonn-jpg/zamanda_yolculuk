import numpy as np
from zt005.amr_interp import (
    prolongate_scalar, prolongate_tensor, prolongate_vector,
    restrict_scalar, inject_boundary,
)


def test_prolongate_constant():
    c = np.ones((5, 5, 5)) * 3.0
    f = prolongate_scalar(c, factor=2)
    assert f.shape == (9, 9, 9)
    assert np.allclose(f, 3.0)


def test_prolongate_linear():
    """Linear field: interpolation should be exact."""
    N = 5
    x = np.linspace(0, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    c = 2.0 * X + 3.0 * Y + 5.0 * Z
    f = prolongate_scalar(c, factor=2)
    Nf = 9
    xf = np.linspace(0, 1, Nf)
    Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
    expected = 2.0 * Xf + 3.0 * Yf + 5.0 * Zf
    assert np.allclose(f, expected, atol=1e-12)


def test_prolongate_tensor_shapes():
    c = np.zeros((4, 4, 4, 3, 3))
    c[..., 0, 0] = 1.0
    f = prolongate_tensor(c, factor=2)
    assert f.shape == (7, 7, 7, 3, 3)
    assert np.allclose(f[..., 0, 0], 1.0)


def test_prolongate_vector_shapes():
    c = np.zeros((4, 4, 4, 3))
    c[..., 0] = 2.0
    f = prolongate_vector(c, factor=2)
    assert f.shape == (7, 7, 7, 3)
    assert np.allclose(f[..., 0], 2.0)


def test_restrict_identity_at_grid():
    """Restriction of prolongated constant is identity."""
    c = np.ones((5, 5, 5)) * 7.0
    f = prolongate_scalar(c, factor=2)
    c2 = restrict_scalar(f, factor=2)
    assert np.allclose(c2, 7.0)


def test_inject_boundary_replaces_edges():
    fine = np.ones((9, 9, 9)) * 10.0
    coarse_prol = np.ones((9, 9, 9)) * 20.0
    out = inject_boundary(fine, coarse_prol, n_layers=2)
    # Interior untouched
    assert out[4, 4, 4] == 10.0
    # Boundary replaced
    assert out[0, 4, 4] == 20.0
    assert out[-1, 4, 4] == 20.0




def test_cubic_prolongate_quadratic():
    """Cubic interpolation exact for quadratic (order 3 polynomial)."""
    N = 7
    x = np.linspace(0, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    c = 3.0 * X**2 + 2.0 * Y**2 + 1.0 * Z**2
    from zt005.amr_interp import prolongate_scalar_cubic
    f = prolongate_scalar_cubic(c, factor=2)
    Nf = (N - 1) * 2 + 1
    xf = np.linspace(0, 1, Nf)
    Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
    expected = 3.0 * Xf**2 + 2.0 * Yf**2 + 1.0 * Zf**2
    # Interior only (boundary clamped stencil degrades)
    m = slice(2, Nf-2)
    assert np.allclose(f[m, m, m], expected[m, m, m], atol=1e-10)
