import numpy as np
from zt005.bssn_dissipation import (
    ko_dissipation, sommerfeld_boundary, apply_dissipation_to_state,
)


def test_ko_constant_zero():
    f = np.ones((10, 10, 10))
    d = ko_dissipation(f, h=0.1)
    # interior should be zero for constant field
    assert np.allclose(d[3:-3, 3:-3, 3:-3], 0.0, atol=1e-14)


def test_ko_linear_interior_zero():
    x = np.arange(10)
    f = np.broadcast_to(x[None, None, :], (10, 10, 10)).astype(float)
    d = ko_dissipation(f, h=0.1)
    assert np.allclose(d[3:-3, 3:-3, 3:-3], 0.0, atol=1e-12)


def test_sommerfeld_outflow_alpha_in_range():
    f = np.zeros((5, 5, 5))
    f[-1, :, :] = 1.0
    f_old = np.zeros_like(f)
    out = sommerfeld_boundary(f, f_old, h=0.1, dt=0.05)
    # should be finite and not equal to old
    assert np.all(np.isfinite(out))
    assert not np.allclose(out, f)


def test_apply_dissipation_preserves_shape():
    state = {
        "phi": np.zeros((5,5,5)),
        "gtilde": np.zeros((5,5,5,3,3)),
        "K": np.zeros((5,5,5)),
        "Atilde": np.zeros((5,5,5,3,3)),
    }
    out = apply_dissipation_to_state(state, h=0.1, dt=0.01)
    for k in state:
        assert out[k].shape == state[k].shape
        assert np.allclose(out[k], 0.0)
