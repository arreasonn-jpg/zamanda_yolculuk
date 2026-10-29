"""Long-run stability test for v2.0 BSSN solver.

Note: Full 1000+ step test is too slow for CI (~36 min on N=7).
This test runs a quick sanity check (100 steps on N=5).
For full validation, run manually with larger N and n_steps.
"""
from __future__ import annotations
import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_v2 import rk4_step_v2
from zt005.bssn import BSSN, hamiltonian_constraint, momentum_constraint


def run_longrun_v2(N=5, L=4.0, A=1e-4, sigma=0.5, dt=1e-4, n_steps=100):
    """Run v2.0 BSSN for n_steps and track constraint norms."""
    s = minkowski_state(N)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2
    s.phi = A * np.exp(-R2 / (2.0 * sigma**2))
    
    h = L / (N - 1)
    
    constraint_history = []
    max_phi_history = []
    max_K_history = []
    max_alpha_history = []
    max_beta_history = []
    
    for step in range(n_steps):
        s = rk4_step_v2(s, h, dt)
        
        # Compute constraints every 20 steps
        if step % 20 == 0 or step == n_steps - 1:
            H_vals = []
            M_vals = []
            for ix in range(N):
                for iy in range(N):
                    for iz in range(N):
                        b = BSSN(
                            float(s.phi[ix,iy,iz]), 
                            s.gtilde[ix,iy,iz],
                            float(s.K[ix,iy,iz]), 
                            s.Atilde[ix,iy,iz],
                            np.zeros(3)
                        )
                        H_vals.append(abs(hamiltonian_constraint(b, 0.0)))
                        M_vals.append(np.linalg.norm(momentum_constraint(b, np.zeros(3))))
            
            constraint_history.append({
                "step": step,
                "mean_H": float(np.mean(H_vals)),
                "max_H": float(np.max(H_vals)),
                "mean_M": float(np.mean(M_vals)),
                "max_M": float(np.max(M_vals)),
            })
            
            max_phi_history.append(float(np.max(np.abs(s.phi))))
            max_K_history.append(float(np.max(np.abs(s.K))))
            max_alpha_history.append(float(np.max(np.abs(s.alpha - 1.0))))
            max_beta_history.append(float(np.max(np.abs(s.beta))))
            
            if step % 20 == 0:
                print(f"Step {step}: max|H|={constraint_history[-1]['max_H']:.2e}, "
                      f"max|phi|={max_phi_history[-1]:.2e}, "
                      f"max|K|={max_K_history[-1]:.2e}")
    
    return {
        "constraint_history": constraint_history,
        "max_phi_history": max_phi_history,
        "max_K_history": max_K_history,
        "max_alpha_history": max_alpha_history,
        "max_beta_history": max_beta_history,
        "final_state": s,
        "params": {"N": N, "L": L, "A": A, "sigma": sigma, "dt": dt, "n_steps": n_steps}
    }


def test_v2_longrun_stability():
    """Quick sanity check: v2.0 BSSN remains stable for 100 steps on N=5."""
    result = run_longrun_v2(N=5, n_steps=100, dt=1e-4)
    
    # Check that constraints remain bounded
    final_H = result["constraint_history"][-1]["max_H"]
    final_M = result["constraint_history"][-1]["max_M"]
    
    print(f"\nFinal max|H| = {final_H:.2e}")
    print(f"Final max|M| = {final_M:.2e}")
    print(f"Final max|phi| = {result['max_phi_history'][-1]:.2e}")
    print(f"Final max|K| = {result['max_K_history'][-1]:.2e}")
    print(f"Final max|alpha-1| = {result['max_alpha_history'][-1]:.2e}")
    print(f"Final max|beta| = {result['max_beta_history'][-1]:.2e}")
    
    # Constraints should stay small (not blow up)
    assert final_H < 1.0, f"Hamiltonian constraint blew up: {final_H}"
    assert final_M < 1.0, f"Momentum constraint blew up: {final_M}"
    
    # Fields should remain finite
    assert np.all(np.isfinite(result["final_state"].phi))
    assert np.all(np.isfinite(result["final_state"].K))
    assert np.all(np.isfinite(result["final_state"].alpha))
    assert np.all(np.isfinite(result["final_state"].beta))
    assert np.all(np.isfinite(result["final_state"].Gamma))
    assert np.all(np.isfinite(result["final_state"].gtilde))
    assert np.all(np.isfinite(result["final_state"].Atilde))
    
    # Alpha should stay positive
    assert np.all(result["final_state"].alpha > 0)
    
    print("\n✅ v2.0 long-run stability test PASSED (100 steps, N=5)")


def test_v2_longrun_manual():
    """Manual test for 1000+ steps - run separately due to time.
    
    Usage: python -m pytest tests/test_bssn_v2_longrun.py::test_v2_longrun_manual -v -s
    """
    import os
    if os.environ.get("RUN_FULL_LONGRUN") != "1":
        import pytest
        pytest.skip("Set RUN_FULL_LONGRUN=1 to run full 1000-step test")
    
    result = run_longrun_v2(N=7, n_steps=1000, dt=2e-4)
    
    final_H = result["constraint_history"][-1]["max_H"]
    final_M = result["constraint_history"][-1]["max_M"]
    
    print(f"\n1000 steps - Final max|H| = {final_H:.2e}")
    print(f"1000 steps - Final max|M| = {final_M:.2e}")
    
    assert final_H < 1.0
    assert final_M < 1.0
    assert np.all(np.isfinite(result["final_state"].phi))
    
    print("\n✅ v2.0 long-run stability test PASSED (1000 steps)")


if __name__ == "__main__":
    test_v2_longrun_stability()
    # test_v2_longrun_manual()  # Run manually with RUN_FULL_LONGRUN=1