"""Unit tests for the v14.0 inspiral tracking helpers."""
import numpy as np

from zt005.amr.amr_bh_track import alpha_pair, find_bh_centers, \
    peak_centroids

N = 17
L = 20.0
x = np.linspace(-L / 2, L / 2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")


class _FakeState:
    """Minimal duck-type: a BSSN state with an Atilde field."""

    def __init__(self, Atilde):
        self.Atilde = Atilde


def _two_peaks(p_left=(-2.0, 0.25), p_right=(2.0, -0.25), w=0.6):
    r_l2 = (X - p_left[0]) ** 2 + (Y - p_left[1]) ** 2 + Z ** 2
    r_r2 = (X - p_right[0]) ** 2 + (Y - p_right[1]) ** 2 + Z ** 2
    A = (np.exp(-r_l2 / (2 * w ** 2)) * 0.3
         + np.exp(-r_r2 / (2 * w ** 2)) * 0.3)
    A = A[..., None, None] * np.ones((N, N, N, 3, 3))
    return _FakeState(A)


def test_peak_centroids_recovers_off_axis_pair():
    s = _two_peaks()
    c_l, c_r, A2 = peak_centroids(s, L, frac=0.05)
    assert A2.shape == (N, N, N)
    # x positions recovered near +/-2
    assert abs(c_l[0] + 2.0) < 0.5
    assert abs(c_r[0] - 2.0) < 0.5
    # y positions recovered near +/-0.25 on opposite signs
    assert abs(c_l[1] - 0.25) < 0.75
    assert abs(c_r[1] + 0.25) < 0.75
    # left is the x<0 one
    assert c_l[0] < 0.0 < c_r[0]


def test_peak_centroids_uses_weighted_mean():
    # strongly asymmetric peak amplitude split should still center
    # each half's centroid inside its own half-space
    r_l2 = (X + 2.0) ** 2 + Y ** 2 + Z ** 2
    r_r2 = (X - 3.0) ** 2 + Y ** 2 + Z ** 2
    A = (np.exp(-r_l2 / (2 * 0.6 ** 2)) * 1.0
         + np.exp(-r_r2 / (2 * 0.6 ** 2)) * 0.1)
    s = _FakeState(A[..., None, None]
                   * np.ones((N, N, N, 3, 3)))
    c_l, c_r, _ = peak_centroids(s, L)
    assert abs(c_l[0] + 2.0) < 0.5
    assert abs(c_r[0] - 3.0) < 0.6


def _two_wells(p_left=(-2.0, 0.25), p_right=(2.0, -0.25), w=1.0):
    r_l2 = (X - p_left[0]) ** 2 + (Y - p_left[1]) ** 2 + Z ** 2
    r_r2 = (X - p_right[0]) ** 2 + (Y - p_right[1]) ** 2 + Z ** 2
    al = 0.9 - 0.3 * np.exp(-r_l2 / (2 * w ** 2)) \
        - 0.3 * np.exp(-r_r2 / (2 * w ** 2))
    s = _FakeState.__new__(_FakeState)
    s.alpha = al
    return s


def test_alpha_pair_finds_half_wells():
    s = _two_wells()
    c_l, c_r, depths = alpha_pair(s, L)
    assert abs(c_l[0] + 2.0) < 0.75
    assert abs(c_r[0] - 2.0) < 0.75
    # well centre sits between coarse nodes; depth is node-quantized
    assert 0.6 <= depths[0] < 0.8
    assert abs(depths[1] - depths[0]) < 1e-9
