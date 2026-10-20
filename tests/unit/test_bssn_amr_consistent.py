"""Constraint-consistent AMR injection (v12.15) unit tests.

rederive_Gamma must reproduce the BSSN identity Gamma^i = gtilde^jk
G^k_ij = -d_j gtilde^ij (exact in the continuum, O(h^2) discrete) for
a smooth metric perturbation gtilde = exp(S) with tr S = 0 (so
det gtilde = 1 exactly).
"""
import numpy as np

from zt005.bssn_grid import d_central
from zt005.amr.bssn_amr_consistent import rederive_Gamma, gamma_consistent_band


class _FakeState:
    """Minimal attribute container for band tests."""

    def __init__(self, Gamma, gtilde=None):
        self.Gamma = Gamma
        if gtilde is None:
            N = Gamma.shape[0]
            gtilde = np.broadcast_to(np.eye(3),
                                     (N, N, N, 3, 3)).copy()
        self.gtilde = gtilde


def _traceless_sym(N):
    x = np.linspace(-1.0, 1.0, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    eps = 1e-3
    S = np.zeros((N, N, N, 3, 3))
    S[..., 0, 0] = eps * (2 * X ** 2 - Y ** 2 - Z ** 2) / 3.0
    S[..., 1, 1] = eps * (2 * Y ** 2 - Z ** 2 - X ** 2) / 3.0
    S[..., 2, 2] = eps * (2 * Z ** 2 - X ** 2 - Y ** 2) / 3.0
    S[..., 0, 1] = S[..., 1, 0] = eps * X * Y
    S[..., 0, 2] = S[..., 2, 0] = eps * X * Z
    S[..., 1, 2] = S[..., 2, 1] = eps * Y * Z
    return S


def test_rederive_matches_neg_gradient_identity():
    """Gamma^i = -d_j gtilde^ij (up to O(h^2)) for a det-1 gtilde."""
    N = 21
    S = _traceless_sym(N)
    # gtilde = exp(S) via eigen-decomposition (S symmetric)
    w, V = np.linalg.eigh(S)
    gtilde = np.einsum('...ij,...j,...kj->...ik', V, np.exp(w), V)
    assert np.allclose(np.linalg.det(gtilde), 1.0, atol=1e-12)

    h = 0.1
    Gam = rederive_Gamma(gtilde, h)

    # Reference: -d_j gtilde^{ij}, same central operator
    gtinv = np.linalg.inv(gtilde)
    ref = np.zeros_like(Gam)
    for i in range(3):
        for j in range(3):
            dx, dy, dz = d_central(gtinv[..., i, j], h)
            ref[..., i] -= np.stack([dx, dy, dz], axis=-1)[..., j]
    scale = max(np.max(np.abs(ref)), 1e-30)
    rel = np.max(np.abs(Gam - ref)) / scale
    assert rel < 0.05, f"identity violated: rel={rel:.4f}"


def test_rederive_zero_for_flat_gtilde():
    N = 8
    gt = np.broadcast_to(np.eye(3), (N, N, N, 3, 3)).copy()
    Gam = rederive_Gamma(gt, 0.5)
    assert np.max(np.abs(Gam)) == 0.0


def test_band_overwrite_keeps_interior():
    N = 10
    old = _FakeState(np.ones((N, N, N, 3)))
    fresh = _FakeState(np.zeros((N, N, N, 3)))
    out = gamma_consistent_band(old, fresh, n_layers=3, h_f=0.5)
    assert out is not old.Gamma
    # interior (middle cells) unchanged
    assert np.all(np.isclose(out[3:N - 3, 3:N - 3, 3:N - 3], 1.0))
    # band cells replaced by rederived Gamma of flat gtilde = 0
    assert np.max(np.abs(out[0, :, :])) == 0.0
    assert np.max(np.abs(out[:, N - 1, :])) == 0.0
