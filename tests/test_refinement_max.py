import numpy as np
import pytest
from zt005.amr.refinement import restrict, restrict_max


def test_restrict_max_constant_matches_restrict():
    """Constant field: both give the same result."""
    N_f = 5
    f = np.ones((N_f, N_f, N_f)) * 3.0
    a = restrict(f, 11, (3, 3, 3))
    b = restrict_max(f, 11, (3, 3, 3))
    valid = ~np.isnan(a)
    assert np.allclose(a[valid], b[valid])


def test_restrict_max_preserves_peak():
    """Sharp peak in center: restrict_max should preserve it better."""
    N_f = 5
    f = np.zeros((N_f, N_f, N_f))
    f[2, 2, 2] = 1.0
    a = restrict(f, 11, (3, 3, 3))
    b = restrict_max(f, 11, (3, 3, 3), blend=1.0)
    # Average: peak / 8
    # Max: peak itself
    assert a[4, 4, 4] == pytest.approx(1.0 / 8.0)
    assert b[4, 4, 4] == pytest.approx(1.0)


def test_restrict_max_blend_default():
    """blend=0.5: halfway between avg and max."""
    N_f = 5
    f = np.zeros((N_f, N_f, N_f))
    f[2, 2, 2] = 1.0
    b = restrict_max(f, 11, (3, 3, 3), blend=0.5)
    assert b[4, 4, 4] == pytest.approx(0.5 * (1.0 / 8.0) + 0.5 * 1.0)


def test_restrict_max_nan_outside():
    f = np.ones((5, 5, 5))
    b = restrict_max(f, 11, (3, 3, 3))
    assert np.isnan(b[0, 0, 0])
    assert not np.isnan(b[3, 3, 3])


def test_restrict_max_shapes():
    f = np.ones((5, 5, 5))
    b = restrict_max(f, 21, (8, 8, 8))
    assert b.shape == (21, 21, 21)
