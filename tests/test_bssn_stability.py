import numpy as np
from zt005.bssn_stability import gaussian_phi_grid, run_stability


def test_gaussian_symmetric():
    g = gaussian_phi_grid(N=11, L=4.0, A=1.0, sigma=0.5)
    # centered -> value at center is max
    c = g[5, 5, 5]
    assert c > g[0, 0, 0]
    assert c > g[-1, -1, -1]


def test_stability_runs():
    out = run_stability(N=7, L=4.0, A=1e-4, sigma=0.5, dt=1e-3, n_steps=3)
    assert out["phi"].shape == (7, 7, 7)
    assert out["constraint_norm"].shape == (7**3,)
    assert np.all(np.isfinite(out["constraint_norm"]))


def test_stability_small_perturbation_bounded():
    out = run_stability(N=7, L=4.0, A=1e-4, sigma=0.5, dt=1e-4, n_steps=3)
    # Hamiltonian constraint stays small for A=1e-4
    assert np.max(out["constraint_norm"]) < 1e-2


def test_stability_max_phi_bounded():
    out = run_stability(N=7, L=4.0, A=1e-4, sigma=0.5, dt=1e-4, n_steps=3)
    assert out["max_phi"] < 1e-3
