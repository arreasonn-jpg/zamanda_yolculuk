#!/usr/bin/env python3
"""Quick test of AMR moving patch with fewer steps."""
import numpy as np
import json
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from zt005.bssn_ic_binary import brill_lindquist_state as brill_lindquist_ic
from zt005.bssn_state import BSSNState
from zt005.bssn_constraints import hamiltonian_constraint_grid
from zt005.amr.amr_projection import project_phi_band_covariant, momentum_band_correct, _band_mask
from zt005.bssn_bc import sommerfeld_all
from zt005.bssn_rk4_beta3 import rk4_step_beta3

def find_bh_centers_alpha_pair(state, dx):
    """Find BH centers using alpha_pair (lapse funnel) tracking - simplified."""
    N = state.phi.shape[0]
    alpha = state.alpha
    
    # Simple approach: find two deepest minima
    # Use a simple peak finding on -alpha
    from scipy.ndimage import minimum_filter
    
    alpha_smooth = alpha.copy()
    footprint = np.ones((3,3,3))
    local_min = (alpha_smooth == minimum_filter(alpha_smooth, footprint=footprint))
    
    min_coords = np.argwhere(local_min)
    min_values = alpha_smooth[local_min]
    
    sorted_idx = np.argsort(min_values)
    min_coords = min_coords[sorted_idx]
    min_values = min_values[sorted_idx]
    
    centers = []
    for coord in min_coords:
        if len(centers) >= 2:
            break
        ok = True
        for c in centers:
            if np.linalg.norm(coord - c) < 5:
                ok = False
                break
        if ok:
            centers.append(coord)
    
    if len(centers) < 2:
        center = N // 2
        sep_cells = int(15.44 / dx)
        centers = [np.array([center - sep_cells//2, center, center]),
                   np.array([center + sep_cells//2, center, center])]
    
    return np.array(centers)

def compute_com(centers, masses=None):
    if masses is None:
        masses = np.ones(len(centers))
    return np.average(centers, axis=0, weights=masses)

def run_moving_patch_test():
    # Configuration - smaller for quick test
    N_c = 33
    L_c = 160.0
    N_f = 25
    L_f = 40.0
    dt = 0.03125
    steps = 10  # Quick test with 10 steps
    kappa = 2.0
    n_layers = 4
    mask_threshold = 0.5
    
    dx_c = L_c / N_c
    dx_f = L_f / N_f
    
    print(f"Configuration:")
    print(f"  Coarse: N={N_c}, L={L_c}, dx={dx_c:.4f}")
    print(f"  Fine:   N={N_f}, L={L_f}, dx={dx_f:.4f}")
    print(f"  dt={dt}, steps={steps}, kappa={kappa}")
    print(f"  Projection: covariant phi + momentum, n_layers={n_layers}, mask_threshold={mask_threshold}")
    
    # Initial data: Brill-Lindquist binary
    m1 = m2 = 0.5
    d = 15.44
    
    print("\nGenerating coarse grid initial data...")
    state_c = brill_lindquist_ic(N_c, L_c, m1, m2, d)
    
    print("Generating fine grid initial data...")
    state_f = brill_lindquist_ic(N_f, L_f, m1, m2, d)
    
    # Initial metrics
    H_c = hamiltonian_constraint_grid(state_c.phi, state_c.gtilde, state_c.K, state_c.Atilde, h=dx_c)
    H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h=dx_f)
    
    band_mask = _band_mask(state_f.phi, dx_f, n_layers, mask_threshold)
    H_f_band = np.max(np.abs(H_f[band_mask])) if np.any(band_mask) else 0.0
    
    centers = find_bh_centers_alpha_pair(state_c, dx_c)
    com = compute_com(centers)
    separation = np.linalg.norm(centers[1] - centers[0]) * dx_c
    
    print(f"\nStep 0:")
    print(f"  H_c_max = {np.max(np.abs(H_c)):.3e}")
    print(f"  H_f_max = {np.max(np.abs(H_f)):.3e}")
    print(f"  H_f_band = {H_f_band:.3e}")
    print(f"  Separation = {separation:.3f}")
    print(f"  COM = {com * dx_c}")
    
    # Time evolution
    for step in range(1, steps + 1):
        # Evolve coarse grid
        old_c = BSSNState(
            phi=state_c.phi.copy(), gtilde=state_c.gtilde.copy(),
            K=state_c.K.copy(), Atilde=state_c.Atilde.copy(),
            Gamma=state_c.Gamma.copy(), alpha=state_c.alpha.copy(),
            beta=state_c.beta.copy()
        )
        state_c = rk4_step_beta3(state_c, dx_c, dt)
        bc_c = sommerfeld_all(state_c, old_c, dx_c, dt)
        state_c = BSSNState(
            phi=bc_c["phi"], gtilde=bc_c["gtilde"],
            K=bc_c["K"], Atilde=bc_c["Atilde"],
            Gamma=bc_c["Gamma"], alpha=bc_c["alpha"],
            beta=bc_c["beta"]
        )
        
        # Prolongate coarse -> fine (recenter fine grid on COM)
        centers = find_bh_centers_alpha_pair(state_c, dx_c)
        com = compute_com(centers)
        
        # Simple prolongation: reinitialize fine grid at new position
        new_d = separation
        state_f = brill_lindquist_ic(N_f, L_f, m1, m2, new_d)
        
        # Evolve fine grid
        old_f = BSSNState(
            phi=state_f.phi.copy(), gtilde=state_f.gtilde.copy(),
            K=state_f.K.copy(), Atilde=state_f.Atilde.copy(),
            Gamma=state_f.Gamma.copy(), alpha=state_f.alpha.copy(),
            beta=state_f.beta.copy()
        )
        state_f = rk4_step_beta3(state_f, dx_f, dt)
        bc_f = sommerfeld_all(state_f, old_f, dx_f, dt)
        state_f = BSSNState(
            phi=bc_f["phi"], gtilde=bc_f["gtilde"],
            K=bc_f["K"], Atilde=bc_f["Atilde"],
            Gamma=bc_f["Gamma"], alpha=bc_f["alpha"],
            beta=bc_f["beta"]
        )
        
        # Projection: covariant phi + momentum correction in band
        H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h=dx_f)
        state_f_phi, _ = project_phi_band_covariant(state_f.phi, state_f.gtilde, dx_f, H_f, 
                                                  n_layers=n_layers, mask_threshold=mask_threshold)
        state_f = BSSNState(
            phi=state_f_phi, gtilde=state_f.gtilde,
            K=state_f.K, Atilde=state_f.Atilde,
            Gamma=state_f.Gamma, alpha=state_f.alpha,
            beta=state_f.beta
        )
        At_new, _, _ = momentum_band_correct(state_f.Atilde, state_f.gtilde, dx_f, 
                                              n_layers=n_layers, passes=1)
        state_f = BSSNState(
            phi=state_f.phi, gtilde=state_f.gtilde,
            K=state_f.K, Atilde=At_new,
            Gamma=state_f.Gamma, alpha=state_f.alpha,
            beta=state_f.beta
        )
        
        # Compute metrics
        H_c = hamiltonian_constraint_grid(state_c.phi, state_c.gtilde, state_c.K, state_c.Atilde, h=dx_c)
        H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h=dx_f)
        H_f_band = np.max(np.abs(H_f[band_mask])) if np.any(band_mask) else 0.0
        
        # Track BH positions
        centers = find_bh_centers_alpha_pair(state_c, dx_c)
        com = compute_com(centers)
        separation = np.linalg.norm(centers[1] - centers[0]) * dx_c
        
        H_c_max = np.max(np.abs(H_c))
        H_f_max = np.max(np.abs(H_f))
        
        print(f"\nStep {step}:")
        print(f"  H_c_max = {H_c_max:.3e}")
        print(f"  H_f_max = {H_f_max:.3e}")
        print(f"  H_f_band = {H_f_band:.3e}")
        print(f"  Separation = {separation:.3f}")
        print(f"  COM = ({com[0]*dx_c:.3f}, {com[1]*dx_c:.3f}, {com[2]*dx_c:.3f})")
        
        # Check for blowup
        if H_c_max > 1e3 or np.isnan(H_c_max):
            print(f"\n*** COARSE GRID BLOWUP at step {step} ***")
            break
        if H_f_max > 1e3 or np.isnan(H_f_max):
            print(f"\n*** FINE GRID BLOWUP at step {step} ***")
            break
    
    print(f"\n=== QUICK TEST SUMMARY ===")
    print(f"Steps completed: {step}")
    print(f"H_c_max: {np.max(np.abs(H_c)):.3e}")
    print(f"H_f_band: {H_f_band:.3e}")
    print(f"Separation: {separation:.3f}")
    
    return True

if __name__ == '__main__':
    run_moving_patch_test()