import numpy as np
from zt005.lqc_scalar_evolution import (
    H_squared_lqc, rho_scalar, evolve_lqc_scalar, RHO_C,
)


def test_h_squared_sign():
    """H^2 >= 0 for rho <= rho_c, < 0 for rho > rho_c."""
    assert H_squared_lqc(0.5 * RHO_C) > 0
    assert H_squared_lqc(RHO_C) == 0.0
    assert H_squared_lqc(1.5 * RHO_C) < 0


def test_rho_scalar_kinetic_only():
    V = lambda phi: 0.0
    assert rho_scalar(1.0, 2.0, V) == 2.0


def test_evolve_static_at_zero():
    """phi=0, phidot=0, V=0 -> nothing moves."""
    V = lambda phi: 0.0
    dV = lambda phi: 0.0
    out = evolve_lqc_scalar(0.0, 0.0, V, dV, dt=1e-3, n_steps=10)
    assert np.allclose(out["phi"], 0.0)
    assert np.allclose(out["rho"], 0.0)


def test_evolve_de_sitter():
    """Constant V>0 -> accelerated expansion (H > 0)."""
    V0 = 1e-10  # small positive potential
    V = lambda phi: V0
    dV = lambda phi: 0.0
    out = evolve_lqc_scalar(0.0, 0.0, V, dV, dt=1e-3, n_steps=100)
    assert out["H"][0] > 0
    assert out["a"][-1] > out["a"][0]


def test_evolve_bounce_near_rho_c():
    """Start near rho_c: Hubble should be small."""
    V0 = RHO_C * 0.49  # 49% of critical
    V = lambda phi: V0
    dV = lambda phi: 0.0
    out = evolve_lqc_scalar(0.0, 0.0, V, dV, dt=1e-3, n_steps=5)
    # H ~ sqrt((8 pi G/3) rho (1 - 0.49)) ~ still positive
    assert out["H"][0] > 0


def test_shapes():
    V = lambda phi: 1e-10
    dV = lambda phi: 0.0
    out = evolve_lqc_scalar(0.0, 0.0, V, dV, dt=1e-3, n_steps=20)
    for k in ("t", "phi", "rho", "H", "a"):
        assert k in out
        assert len(out[k]) == 21
