"""Schwarzschild puncture evolution test with v2.0-beta3 (shift + BC)."""
from __future__ import annotations
import numpy as np
from zt005.bssn_ic_puncture import schwarzschild_puncture_state
from zt005.bssn_state import BSSNState
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full


def test_schwarzschild_puncture_beta3_evolution():
    """Evolve Schwarzschild puncture with beta3 driver for 50 steps (CI-friendly)."""
    N = 9
    L = 16.0
    h = L / (N - 1)
    dt = 2e-4
    n_steps = 50
    
    # Create puncture initial data
    s_base = schwarzschild_puncture_state(N, L, M=1.0, eps=1.0)
    
    # Convert to BSSNState (beta3 expects BSSNState)
    s = BSSNState(
        phi=s_base.phi.copy(),
        gtilde=s_base.gtilde.copy(),
        K=s_base.K.copy(),
        Atilde=s_base.Atilde.copy(),
        Gamma=s_base.Gamma.copy(),
        alpha=s_base.alpha.copy(),
        beta=s_base.beta.copy(),
    )
    
    # Track constraint
    H_history = []
    
    for step in range(n_steps):
        s = rk4_step_beta3(s, h, dt)
        
        if step % 20 == 0 or step == n_steps - 1:
            H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
            # Exterior region only (r > 5*eps)
            x = np.linspace(-L/2, L/2, N)
            X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
            R = np.sqrt(X**2 + Y**2 + Z**2)
            mask = R > 5.0
            H_ext = float(np.max(np.abs(H[mask])))
            H_history.append(H_ext)
            print(f"Step {step}: max|H|_ext = {H_ext:.2e}")
    
    # Constraint should not grow significantly
    growth = H_history[-1] / H_history[0] if H_history[0] > 0 else 1.0
    print(f"\nH growth factor: {growth:.4f}")
    
    # Should remain bounded (growth < 1.1 for 100 steps)
    assert growth < 1.1, f"Hamiltonian constraint grew by factor {growth}"
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(s.K))
    assert np.all(np.isfinite(s.alpha))
    assert np.all(np.isfinite(s.beta))
    assert np.all(s.alpha > 0)
    
    print("\n✅ Schwarzschild puncture evolution with beta3 PASSED")


def test_schwarzschild_puncture_beta3_longrun():
    """Longer run (500 steps) - run manually."""
    import os
    if os.environ.get("RUN_FULL_LONGRUN") != "1":
        import pytest
        pytest.skip("Set RUN_FULL_LONGRUN=1 to run full 500-step test")
    
    N = 11
    L = 20.0
    h = L / (N - 1)
    dt = 5e-5
    n_steps = 500
    
    s_base = schwarzschild_puncture_state(N, L, M=1.0, eps=1.0)
    s = BSSNState(
        phi=s_base.phi.copy(),
        gtilde=s_base.gtilde.copy(),
        K=s_base.K.copy(),
        Atilde=s_base.Atilde.copy(),
        Gamma=s_base.Gamma.copy(),
        alpha=s_base.alpha.copy(),
        beta=s_base.beta.copy(),
    )
    
    H_history = []
    
    for step in range(n_steps):
        s = rk4_step_beta3(s, h, dt)
        
        if step % 100 == 0 or step == n_steps - 1:
            H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
            x = np.linspace(-L/2, L/2, N)
            X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
            R = np.sqrt(X**2 + Y**2 + Z**2)
            mask = R > 5.0
            H_ext = float(np.max(np.abs(H[mask])))
            H_history.append(H_ext)
            print(f"Step {step}: max|H|_ext = {H_ext:.2e}")
    
    growth = H_history[-1] / H_history[0] if H_history[0] > 0 else 1.0
    print(f"\n500-step H growth factor: {growth:.4f}")
    
    assert growth < 1.05, f"Hamiltonian constraint grew by factor {growth}"
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(s.K))
    assert np.all(np.isfinite(s.alpha))
    assert np.all(np.isfinite(s.beta))
    assert np.all(s.alpha > 0)
    
    print("\n✅ Schwarzschild puncture 500-step evolution with beta3 PASSED")


if __name__ == "__main__":
    test_schwarzschild_puncture_beta3_evolution()
    # test_schwarzschild_puncture_beta3_longrun()  # Run manually with RUN_FULL_LONGRUN=1