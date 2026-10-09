import numpy as np
from zt005.bssn_rk4 import evolve_rk4
from zt005.bssn_grid import christoffel_from_grid, dG_from_grid


def test_christoffel_grid_flat_zero():
    N = 5
    gt = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    G, _ = christoffel_from_grid(gt, 0.1)
    assert np.allclose(G, 0.0)


def test_dG_grid_flat_zero():
    N = 5
    G = np.zeros((N,N,N,3,3,3))
    dG = dG_from_grid(G, 0.1)
    assert np.allclose(dG, 0.0)


def test_rk4_minkowski_flat():
    N = 3
    phi = np.zeros((N,N,N))
    gt = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    K = np.zeros((N,N,N))
    At = np.zeros((N,N,N,3,3))
    out = evolve_rk4(phi, gt, K, At, h=0.1, dt=0.001, n_steps=3)
    assert np.allclose(out["phi"], 0.0, atol=1e-12)
    assert np.allclose(out["K"], 0.0, atol=1e-12)


def test_rk4_constraint_does_not_blow_up():
    N = 3
    phi = np.zeros((N,N,N))
    gt = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    K = np.ones((N,N,N)) * 1e-6
    At = np.zeros((N,N,N,3,3))
    out = evolve_rk4(phi, gt, K, At, h=0.1, dt=1e-4, n_steps=4)
    n = out["constraint_norm"]
    assert len(n) == 4
    assert n[-1] < 1e-3
