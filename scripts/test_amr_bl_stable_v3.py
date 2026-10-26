#!/usr/bin/env python3
"""
AMR test with Brill-Lindquist initial data (time-symmetric, K=0)
More stable than Bowen-York for long runs.
Uses the correct BSSNState structure from bssn_state.py
"""

import numpy as np
from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.amr_long_run import AMRConfig
from zt005.bssn_state import BSSNState
from zt005.bssn_rhs_vec import bssn_rhs_vec
from zt005.bssn_constraints import hamiltonian_constraint_grid
from zt005.amr.amr_projection import project_phi_band_covariant, momentum_band_correct, _band_mask
import os

# Stable config from previous tests: N_c=33, L_c=80, N_f=25, L_f=40, dt=0.03125, kappa=2.0
# But use Brill-Lindquist initial data instead of Bowen-York

config = AMRConfig(
    N_c=33,
    L_c=80.0,
    N_f=25,
    L_f=40.0,
    dt=0.03125,
    n_steps=30,
    kappa=2.0,
    # Projection settings
    use_projection=True,
    projection_type="covariant",
    n_layers=4,
    mask_threshold=0.5,
    use_momentum_correction=True,
    momentum_passes=1,
    # Output
    output_every=1,
    output_dir="results/amr_bl_stable",
)

print("=" * 60)
print("AMR Test: Brill-Lindquist Binary BH (Time-Symmetric)")
print("=" * 60)
print(f"Config: N_c={config.N_c}, L_c={config.L_c}, h_c={config.L_c/config.N_c:.3f}")
print(f"        N_f={config.N_f}, L_f={config.L_f}, h_f={config.L_f/config.N_f:.3f}")
print(f"        dt={config.dt}, steps={config.n_steps}, kappa={config.kappa}")
print(f"        Projection: phi={config.use_projection}, momentum={config.use_momentum_correction}")
print(f"        n_layers={config.n_layers}, mask_threshold={config.mask_threshold}")
print(f"        projection_type={config.projection_type}, momentum_passes={config.momentum_passes}")
print("=" * 60)

def run_amr_bl(config: AMRConfig):
    """Run AMR with Brill-Lindquist initial data."""
    os.makedirs(config.output_dir, exist_ok=True)
    
    # Create grids (inline like amr_long_run.py)
    x_c = np.linspace(-config.L_c / 2, config.L_c / 2, config.N_c)
    X_c, Y_c, Z_c = np.meshgrid(x_c, x_c, x_c, indexing="ij")
    h_c = config.L_c / (config.N_c - 1)
    
    x_f = np.linspace(-config.L_f / 2, config.L_f / 2, config.N_f)
    X_f, Y_f, Z_f = np.meshgrid(x_f, x_f, x_f, indexing="ij")
    h_f = config.L_f / (config.N_f - 1)
    
    # Brill-Lindquist parameters (time-symmetric, K=0)
    # Two equal mass BHs at +/- d/2 on x-axis
    d = 12.0  # separation
    m1 = 0.5  # mass each
    m2 = 0.5
    
    # Coarse grid initial data using brill_lindquist_state
    state_c = brill_lindquist_state(config.N_c, config.L_c, m1, m2, d)
    
    # Fine grid initial data
    state_f = brill_lindquist_state(config.N_f, config.L_f, m1, m2, d)
    
    # Band mask for projection
    mask = _band_mask(state_f.phi, h_f, config.n_layers, config.mask_threshold)
    band_cells = np.sum(mask)
    print(f"Band mask cells: {band_cells}")
    
    # Initial constraints
    H_c = hamiltonian_constraint_grid(state_c.phi, state_c.gtilde, state_c.K, state_c.Atilde, h_c)
    H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h_f)
    print(f"Initial H_c max: {np.max(np.abs(H_c)):.3e}")
    print(f"Initial H_f max: {np.max(np.abs(H_f)):.3e}")
    print(f"Initial H_f band max: {np.max(np.abs(H_f[mask])):.3e}")
    
    # Time evolution
    for step in range(config.n_steps):
        t = step * config.dt
        
        # Coarse grid evolution
        rhs_c = bssn_rhs_vec(
            state_c.phi, state_c.gtilde, state_c.K, state_c.Atilde,
            state_c.alpha, h_c,
            sponge=None, kappa=config.kappa,
            beta=state_c.beta, Gamma=state_c.Gamma
        )
        state_c.phi += config.dt * rhs_c["dphi"]
        state_c.K += config.dt * rhs_c["dK"]
        state_c.gtilde += config.dt * rhs_c["dgt"]
        state_c.Atilde += config.dt * rhs_c["dAt"]
        state_c.Gamma += config.dt * rhs_c["dGamma"]
        # alpha and beta are not evolved in this simple test
        
        # Fine grid evolution
        rhs_f = bssn_rhs_vec(
            state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde,
            state_f.alpha, h_f,
            sponge=None, kappa=config.kappa,
            beta=state_f.beta, Gamma=state_f.Gamma
        )
        state_f.phi += config.dt * rhs_f["dphi"]
        state_f.K += config.dt * rhs_f["dK"]
        state_f.gtilde += config.dt * rhs_f["dgt"]
        state_f.Atilde += config.dt * rhs_f["dAt"]
        state_f.Gamma += config.dt * rhs_f["dGamma"]
        
        # Projection every step
        if config.use_projection:
            if config.projection_type == "covariant":
                state_f.phi, dphi_max = project_phi_band_covariant(
                    state_f.phi, state_f.gtilde, h_f, H_f,
                    n_layers=config.n_layers,
                    mask_threshold=config.mask_threshold
                )
            else:
                from zt005.amr.amr_projection import project_phi_band
                state_f.phi, dphi_max = project_phi_band(
                    state_f.phi, h_f, H_f,
                    n_layers=config.n_layers,
                    mask_threshold=config.mask_threshold
                )
        
        if config.use_momentum_correction:
            for _ in range(config.momentum_passes):
                state_f.Atilde, l_before, l_after = momentum_band_correct(
                    state_f.Atilde, state_f.gtilde, h_f,
                    n_layers=config.n_layers,
                    passes=1
                )
        
        # Constraints
        H_c = hamiltonian_constraint_grid(state_c.phi, state_c.gtilde, state_c.K, state_c.Atilde, h_c)
        H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h_f)
        H_f_band = np.max(np.abs(H_f[mask])) if band_cells > 0 else 0
        
        # Momentum constraint
        from zt005.amr.amr_projection import atilde_up_divergence
        l_vec = atilde_up_divergence(state_f.Atilde, state_f.gtilde, h_f)
        l_norm = np.sqrt(l_vec[...,0]**2 + l_vec[...,1]**2 + l_vec[...,2]**2)
        l_before = np.max(l_norm[mask]) if band_cells > 0 else 0
        
        if config.use_momentum_correction:
            for _ in range(config.momentum_passes):
                state_f.Atilde, _, _ = momentum_band_correct(
                    state_f.Atilde, state_f.gtilde, h_f,
                    n_layers=config.n_layers,
                    passes=1
                )
            l_vec = atilde_up_divergence(state_f.Atilde, state_f.gtilde, h_f)
            l_norm = np.sqrt(l_vec[...,0]**2 + l_vec[...,1]**2 + l_vec[...,2]**2)
        l_after = np.max(l_norm[mask]) if band_cells > 0 else 0
        
        # Track BH separation (find minima of phi)
        # Simplified: track coordinate separation
        sep = d  # placeholder
        
        if step % config.output_every == 0:
            print(f"Step {step:3d} t={t:5.2f} | H_c_max={np.max(np.abs(H_c)):.3e} H_f_max={np.max(np.abs(H_f)):.3e} H_f_band={H_f_band:.3e} l_before={l_before:.3e} l_after={l_after:.3e}")
        
        # Check for blowup
        if np.max(np.abs(H_c)) > 1.0 or np.max(np.abs(H_f)) > 1.0:
            print(f"BLOWUP at step {step}: H_c={np.max(np.abs(H_c)):.3e}, H_f={np.max(np.abs(H_f)):.3e}")
            break
    
    # Final summary
    H_c = hamiltonian_constraint_grid(state_c.phi, state_c.gtilde, state_c.K, state_c.Atilde, h_c)
    H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h_f)
    H_f_band = np.max(np.abs(H_f[mask])) if band_cells > 0 else 0
    
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Final H_f_band: {H_f_band:.3e}")
    print(f"Final H_c_max:  {np.max(np.abs(H_c)):.3e}")
    print(f"Steps completed: {step + 1}")
    if np.max(np.abs(H_c)) < 1.0:
        print("[OK] COARSE GRID STABLE")
    else:
        print("[FAIL] COARSE GRID BLOWUP")
    
    return step + 1, np.max(np.abs(H_c)), H_f_band

if __name__ == "__main__":
    run_amr_bl(config)