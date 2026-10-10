import numpy as np
from zt005.bssn_rhs_grid import (
    grad_scalar, hess_scalar, grad_tensor, bssn_rhs_grid,
)


def test_grad_linear():
    x = np.arange(11)
    f = np.broadcast_to(x[None,None,:], (3,3,11)).astype(float)
    dx, dy, dz = grad_scalar(f, h=1.0)
    assert np.allclose(dz[1,1,1:-1], 1.0)


def test_hess_quadratic():
    x = np.arange(11)
    f = np.broadcast_to((x**2)[None,None,:], (3,3,11)).astype(float)
    H = hess_scalar(f, h=1.0)
    assert H.shape == (3, 3, 3, 3, 11)


def test_grad_tensor_identity_zero():
    N = 5
    A = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    dA = grad_tensor(A, h=0.1)
    assert np.allclose(dA, 0.0)


def test_rhs_grid_minkowski():
    N = 5
    phi = np.zeros((N,N,N))
    gt = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    K = np.zeros((N,N,N))
    At = np.zeros((N,N,N,3,3))
    rhs = bssn_rhs_grid(phi, gt, K, At, h=0.1)
    for k in ("dphi", "dgt", "dK", "dAt"):
        assert k in rhs
    assert np.allclose(rhs["dphi"], 0.0)
    assert np.allclose(rhs["dgt"], 0.0)
