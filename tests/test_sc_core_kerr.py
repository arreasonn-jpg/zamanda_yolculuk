"""SC regular core + Kerr (rotating SC) evolution test (v3.0-alpha6)."""
from __future__ import annotations
import numpy as np
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_rk4_mp import rk4_step_mp
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.sc_core_kerr import sc_core_kerr_state


def test_sc_core_kerr_evolution():
    """Evolve rotating SC regular core with moving puncture for 100 steps."""
    N = 11
    L = 40.0
    M = 1.0
    rho0 = 1e-3
    a = 0.3  # slow rotation
    h = L / (N - 1)
    dt = 2e-4
    n_steps = 100  # CI-friendly
    
    s, params = sc_core_kerr_state(N, L, M=M, rho0=rho0, a=a)
    
    print(f"SC Core Kerr params: Rs={params['Rs']:.2e}, R0={params['R0']:.2e}, ell={params['ell']:.2e}, a={a}")
    
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    R0 = params["R0"]
    mask = (R > R0) & (R < L/2)
    
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
    
    print("\n[PASS] SC regular core + Kerr (rotating SC) 100-step PASSED")


def test_sc_core_kerr_evolution_longrun():
    """Full 500-step rotating SC core evolution - run manually."""
    import os
    if os.environ.get("RUN_FULL_LONGRUN") != "1":
        import pytest
        pytest.skip("Set RUN_FULL_LONGRUN=1 to run full 500-step test")
    
    N = 11
    L = 40.0
    M = 1.0
    rho0 = 1e-3
    a = 0.3
    h = L / (N - 1)
    dt = 2e-4
    n_steps = 500
    
    s, params = sc_core_kerr_state(N, L, M=M, rho0=rho0, a=a)
    
    print(f"SC Core Kerr params: Rs={params['Rs']:.2e}, R0={params['R0']:.2e}, ell={params['ell']:.2e}, a={a}")
    
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    R0 = params["R0"]
    mask = (R > R0) & (R < L/2)
    
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
    
    print("\n[PASS] SC regular core + Kerr (rotating SC) 500-step PASSED")


if __name__ == "__main__":
    test_sc_core_kerr_evolution()
    # test_sc_core_kerr_evolution_longrun()  # Run manually with RUN_FULL_LONGRUN=1