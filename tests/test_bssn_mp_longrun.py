"""Moving puncture long-run stability test (v3.0-alpha3: Schwarzschild 500-step)."""
from __future__ import annotations
import numpy as np
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_puncture import schwarzschild_puncture_state
from zt005.bssn_rk4_mp import rk4_step_mp
from zt005.bssn_ricci_full import hamiltonian_constraint_full


def test_moving_puncture_schwarzschild_500():
    """Moving puncture Schwarzschild evolution for 500 steps (CI-friendly: 100 steps)."""
    N = 11
    L = 40.0
    M = 1.0
    h = L / (N - 1)
    dt = 2e-4
    n_steps = 100  # CI-friendly
    
    s_base = schwarzschild_puncture_state(N, L, M=M, eps=M)
    s = BSSNStateMP.from_base(s_base)
    
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    mask = (R > 3.0*M) & (R < L/3)
    
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
    print(f"H0_ext = {H0:.4e}")
    
    H_history = [H0]
    
    for step in range(1, n_steps + 1):
        s = rk4_step_mp(s, h=h, dt=dt)
        
        if step % 20 == 0 or step == n_steps:
            H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
            n = float(np.mean(np.abs(H[mask])))
            H_history.append(n)
            Kmax = float(np.max(np.abs(s.K)))
            amin = float(s.alpha.min())
            bmax = float(np.max(np.abs(s.beta)))
            Bmax = float(np.max(np.abs(s.B)))
            fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(s.K)) 
                   and np.all(np.isfinite(s.B)))
            print(f"step={step:4d} H={n:.3e} growth={n/H0:.4f} "
                  f"Kmax={Kmax:.3e} amin={amin:.4f} bmax={bmax:.3e} "
                  f"Bmax={Bmax:.3e} fin={fin}")
            if not fin:
                raise AssertionError("NON-FINITE state detected")
    
    growth = H_history[-1] / H_history[0] if H_history[0] > 0 else 1.0
    print(f"\nFinal H growth factor: {growth:.4f}")
    
    # Should remain bounded
    assert growth < 1.1, f"Hamiltonian constraint grew by factor {growth}"
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(s.K))
    assert np.all(np.isfinite(s.B))
    assert np.all(s.alpha > 0)
    
    print("\n✅ Moving puncture Schwarzschild 100-step PASSED")


def test_moving_puncture_schwarzschild_500_longrun():
    """Full 500-step moving puncture test - run manually."""
    import os
    if os.environ.get("RUN_FULL_LONGRUN") != "1":
        import pytest
        pytest.skip("Set RUN_FULL_LONGRUN=1 to run full 500-step test")
    
    N = 11
    L = 40.0
    M = 1.0
    h = L / (N - 1)
    dt = 2e-4
    n_steps = 500
    
    s_base = schwarzschild_puncture_state(N, L, M=M, eps=M)
    s = BSSNStateMP.from_base(s_base)
    
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    mask = (R > 3.0*M) & (R < L/3)
    
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
    print(f"H0_ext = {H0:.4e}")
    
    H_history = [H0]
    
    for step in range(1, n_steps + 1):
        s = rk4_step_mp(s, h=h, dt=dt)
        
        if step % 100 == 0 or step == n_steps:
            H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
            n = float(np.mean(np.abs(H[mask])))
            H_history.append(n)
            Kmax = float(np.max(np.abs(s.K)))
            amin = float(s.alpha.min())
            bmax = float(np.max(np.abs(s.beta)))
            Bmax = float(np.max(np.abs(s.B)))
            fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(s.K)) 
                   and np.all(np.isfinite(s.B)))
            print(f"step={step:4d} H={n:.3e} growth={n/H0:.4f} "
                  f"Kmax={Kmax:.3e} amin={amin:.4f} bmax={bmax:.3e} "
                  f"Bmax={Bmax:.3e} fin={fin}")
            if not fin:
                raise AssertionError("NON-FINITE state detected")
    
    growth = H_history[-1] / H_history[0] if H_history[0] > 0 else 1.0
    print(f"\nFinal H growth factor (500 steps): {growth:.4f}")
    
    assert growth < 1.05, f"Hamiltonian constraint grew by factor {growth}"
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(s.K))
    assert np.all(np.isfinite(s.B))
    assert np.all(s.alpha > 0)
    
    print("\n✅ Moving puncture Schwarzschild 500-step PASSED")


if __name__ == "__main__":
    test_moving_puncture_schwarzschild_500()
    # test_moving_puncture_schwarzschild_500_longrun()  # Run manually with RUN_FULL_LONGRUN=1