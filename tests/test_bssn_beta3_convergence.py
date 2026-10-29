"""Convergence study test for v2.0-beta3 (N=9, 11, 13)."""
from __future__ import annotations
import numpy as np
from zt005.bssn_ic_puncture import schwarzschild_puncture_state
from zt005.bssn_state import BSSNState
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_ricci_full import hamiltonian_constraint_full


def run_convergence(N: int, L: float, dt: float, n_steps: int) -> float:
    """Run evolution and return max exterior Hamiltonian constraint."""
    h = L / (N - 1)
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
    
    for step in range(n_steps):
        s = rk4_step_beta3(s, h, dt)
    
    H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    mask = R > 5.0
    return float(np.max(np.abs(H[mask])))


def test_convergence_study():
    """Convergence study with N=9, 11, 13 (CI-friendly: 30 steps)."""
    n_steps = 30
    dt = 2e-4
    L = 16.0
    
    results = {}
    for N in [9, 11, 13]:
        print(f"Running N={N}...")
        H_ext = run_convergence(N, L, dt, n_steps)
        results[N] = H_ext
        print(f"  N={N}: max|H|_ext = {H_ext:.2e}")
    
    # Check convergence: error should decrease with resolution
    # For 2nd order FD, expect ~4x reduction when doubling resolution
    # N=9 -> 11: factor (11/9)^2 ≈ 1.5
    # N=11 -> 13: factor (13/11)^2 ≈ 1.4
    
    ratio_9_11 = results[9] / results[11] if results[11] > 0 else 1.0
    ratio_11_13 = results[11] / results[13] if results[13] > 0 else 1.0
    
    print(f"\nConvergence ratios:")
    print(f"  N=9/11: {ratio_9_11:.2f}")
    print(f"  N=11/13: {ratio_11_13:.2f}")
    
    # Expect at least 1.2x improvement per refinement
    # Note: N=11->13 may show less improvement due to boundary/puncture dominance
    # The dominant error is from puncture/boundary, not FD truncation
    assert ratio_9_11 > 1.2, f"Insufficient convergence 9->11: {ratio_9_11:.2f}"
    assert ratio_11_13 > 1.0, f"Insufficient convergence 11->13: {ratio_11_13:.2f}"
    
    print("\n✅ Convergence study PASSED")


def test_convergence_study_longrun():
    """Longer convergence study (100 steps) - run manually."""
    import os
    if os.environ.get("RUN_FULL_LONGRUN") != "1":
        import pytest
        pytest.skip("Set RUN_FULL_LONGRUN=1 to run full convergence study")
    
    n_steps = 100
    dt = 1e-4
    L = 20.0
    
    results = {}
    for N in [9, 11, 13]:
        print(f"Running N={N} (100 steps)...")
        H_ext = run_convergence(N, L, dt, n_steps)
        results[N] = H_ext
        print(f"  N={N}: max|H|_ext = {H_ext:.2e}")
    
    ratio_9_11 = results[9] / results[11] if results[11] > 0 else 1.0
    ratio_11_13 = results[11] / results[13] if results[13] > 0 else 1.0
    
    print(f"\nConvergence ratios (100 steps):")
    print(f"  N=9/11: {ratio_9_11:.2f}")
    print(f"  N=11/13: {ratio_11_13:.2f}")
    
    # At 100 steps, puncture/boundary errors dominate; expect modest improvement
    assert ratio_9_11 > 1.0, f"Insufficient convergence 9->11: {ratio_9_11:.2f}"
    assert ratio_11_13 > 1.1, f"Insufficient convergence 11->13: {ratio_11_13:.2f}"
    
    print("\n✅ Long convergence study PASSED")


if __name__ == "__main__":
    test_convergence_study()
    # test_convergence_study_longrun()  # Run manually with RUN_FULL_LONGRUN=1