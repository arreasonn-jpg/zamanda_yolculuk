import numpy as np
from zt005.bssn_evolve import evolve_bssn
from zt005.bssn_grid import d_central, dgtilde_from_grid


def test_d_central_linear():
    x = np.linspace(0, 1, 11)
    f = 3.0 * x
    fgrid = np.broadcast_to(f[None, None, :], (3, 3, 11)).copy()
    _, _, dz = d_central(fgrid, x[1]-x[0])
    assert np.allclose(dz[1, 1, 1:-1], 3.0, atol=1e-10)


def test_dgtilde_identity_zero():
    N = 5
    gt = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    dg = dgtilde_from_grid(gt, 0.1)
    assert np.allclose(dg, 0.0)


def test_evolve_minkowski_stationary():
    N = 3
    phi0 = np.zeros((N,N,N))
    gt0 = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    K0 = np.zeros((N,N,N))
    At0 = np.zeros((N,N,N,3,3))
    out = evolve_bssn(phi0, gt0, K0, At0, h=0.1, dt=0.01, n_steps=3)
    assert np.allclose(out["phi"], 0.0, atol=1e-12)
    assert np.allclose(out["K"], 0.0, atol=1e-12)
    assert out["constraint_norm"].shape == (3,)
    assert np.all(out["constraint_norm"] < 1e-20)


def test_evolve_constraint_recorded():
    N = 3
    phi0 = np.zeros((N,N,N))
    gt0 = np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy()
    K0 = np.ones((N,N,N)) * 1e-4
    At0 = np.zeros((N,N,N,3,3))
    out = evolve_bssn(phi0, gt0, K0, At0, h=0.1, dt=0.001, n_steps=4)
    assert len(out["constraint_norm"]) == 4
