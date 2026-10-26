#!/usr/bin/env python3
"""
AMR test with Brill-Lindquist initial data (time-symmetric, K=0)
More stable than Bowen-York for long runs.
"""

import numpy as np
from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.amr_long_run import run_amr_long, AMRConfig

# Stable config from previous tests: N_c=33, L_c=80, N_f=25, L_f=40, dt=0.03125, kappa=2.0
# But use Brill-Lindquist initial data instead of Bowen-York

config = AMRConfig(
    N_coarse=33,
    L_coarse=80.0,
    N_fine=25,
    L_fine=40.0,
    dt=0.03125,
    steps=30,
    kappa=2.0,
    # Projection settings
    project_phi=True,
    project_momentum=True,
    n_layers=4,
    mask_threshold=0.5,
    covariant_phi=True,
    momentum_passes=1,
    # Output
    output_every=1,
    output_dir="results/amr_bl_stable",
)

print("=" * 60)
print("AMR Test: Brill-Lindquist Binary BH (Time-Symmetric)")
print("=" * 60)
print(f"Config: N_c={config.N_coarse}, L_c={config.L_coarse}, h_c={config.L_coarse/config.N_coarse:.3f}")
print(f"        N_f={config.N_fine}, L_f={config.L_fine}, h_f={config.L_fine/config.N_fine:.3f}")
print(f"        dt={config.dt}, steps={config.steps}, kappa={config.kappa}")
print(f"        Projection: phi={config.project_phi}, momentum={config.project_momentum}")
print(f"        n_layers={config.n_layers}, mask_threshold={config.mask_threshold}")
print(f"        covariant_phi={config.covariant_phi}, momentum_passes={config.momentum_passes}")
print("=" * 60)

# Run with Brill-Lindquist initial data
# We need to modify the run_amr_long to accept custom initial data
# For now, let's check if there's a way to pass initial data

# Actually, let's create a modified version that uses Brill-Lindquist
from zt005.bssn import BSSNState, BSSNParams
from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.amr import AMRGrid
from zt005.amr_projection import project_phi_band_covariant, momentum_band_correct
from zt005.bssn_rhs import bssn_rhs
from zt005.bssn_constraints import hamiltonian_constraint, momentum_constraint
import os

def run_amr_bl(config: AMRConfig):
    """Run AMR with Brill-Lindquist initial data."""
    os.makedirs(config.output_dir, exist_ok=True)
    
    # Create grids
    coarse = AMRGrid(config.N_coarse, config.L_coarse)
    fine = AMRGrid(config.N_fine, config.L_fine)
    
    # Brill-Lindquist parameters (time-symmetric, K=0)
    # Two equal mass BHs at +/- d/2 on x-axis
    d = 12.0  # separation
    m1 = 0.5  # mass each
    m2 = 0.5
    
    # Coarse grid initial data using brill_lindquist_state
    state_c = brill_lindquist_state(config.N_coarse, config.L_coarse, m1, m2, d)
    
    # Fine grid initial data
    state_f = brill_lindquist_state(config.N_fine, config.L_fine, m1, m2, d)
    
    params = BSSNParams(
        eta=config.kappa,
        chi=1.0,
        Gamma_driver=1.0,
        dissipation=0.1,
    )
    
    # Band mask for projection
    from src.zt005.amr_projection import band_mask
    mask = band_mask(fine, coarse, n_layers=config.n_layers, threshold=config.mask_threshold)
    band_cells = np.sum(mask)
    print(f"Band mask cells: {band_cells}")
    
    # Initial constraints
    H_c = hamiltonian_constraint(state_c, coarse)
    H_f = hamiltonian_constraint(state_f, fine)
    print(f"Initial H_c max: {np.max(np.abs(H_c)):.3e}")
    print(f"Initial H_f max: {np.max(np.abs(H_f)):.3e}")
    print(f"Initial H_f band max: {np.max(np.abs(H_f[mask])):.3e}")
    
    # Time evolution
    for step in range(config.steps):
        t = step * config.dt
        
        # Coarse grid evolution
        rhs_c = bssn_rhs(state_c, coarse, params)
        for attr in ['phi', 'K', 'Atilde_xx', 'Atilde_xy', 'Atilde_xz', 
                     'Atilde_yy', 'Atilde_yz', 'Atilde_zz',
                     'Gamma_x', 'Gamma_y', 'Gamma_z',
                     'beta_x', 'beta_y', 'beta_z',
                     'B_x', 'B_y', 'B_z']:
            val = getattr(state_c, attr)
            val += config.dt * getattr(rhs_c, attr)
        
        # Fine grid evolution
        rhs_f = bssn_rhs(state_f, fine, params)
        for attr in ['phi', 'K', 'Atilde_xx', 'Atilde_xy', 'Atilde_xz', 
                     'Atilde_yy', 'Atilde_yz', 'Atilde_zz',
                     'Gamma_x', 'Gamma_y', 'Gamma_z',
                     'beta_x', 'beta_y', 'beta_z',
                     'B_x', 'B_y', 'B_z']:
            val = getattr(state_f, attr)
            val += config.dt * getattr(rhs_f, attr)
        
        # Projection every step
        if config.project_phi:
            if config.covariant_phi:
                project_phi_band_covariant(state_f, fine, coarse, mask, n_iter=5)
            else:
                from src.zt005.amr_projection import project_phi_band
                project_phi_band(state_f, fine, coarse, mask, n_iter=5)
        
        if config.project_momentum:
            for _ in range(config.momentum_passes):
                momentum_band_correct(state_f, fine, mask)
        
        # Constraints
        H_c = hamiltonian_constraint(state_c, coarse)
        H_f = hamiltonian_constraint(state_f, fine)
        H_f_band = np.max(np.abs(H_f[mask])) if band_cells > 0 else 0
        
        # Momentum constraint
        M_f = momentum_constraint(state_f, fine)
        M_f_norm = np.sqrt(M_f[0]**2 + M_f[1]**2 + M_f[2]**2)
        l_before = np.max(M_f_norm[mask]) if band_cells > 0 else 0
        
        if config.project_momentum:
            for _ in range(config.momentum_passes):
                momentum_band_correct(state_f, fine, mask)
            M_f = momentum_constraint(state_f, fine)
            M_f_norm = np.sqrt(M_f[0]**2 + M_f[1]**2 + M_f[2]**2)
        l_after = np.max(M_f_norm[mask]) if band_cells > 0 else 0
        
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
    H_c = hamiltonian_constraint(state_c, coarse)
    H_f = hamiltonian_constraint(state_f, fine)
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