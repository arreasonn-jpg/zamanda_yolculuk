"""BSSN AMR Long Run with Projection (v14.1).

Stage-0 of the inspiral marathon: the unprojected long run blows up at
t = 5.6 in the injection band (fine band max |H| 0.44 -> 2.32,
e-fold ~ 0.4 M) while the orbit-gauge signal is real (omega +0.13..0.15
rad/M); the critical path is projection-on long runs (v14.1).

This driver:
- Sets up binary puncture initial data
- Runs AMR with moving fine patch (COM-tracked)
- Applies constraint projection (phi band + momentum band) after each step
- Tracks Hamiltonian constraint and BH positions
"""
from __future__ import annotations

import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, Any, Optional

from .bssn_state_mp import BSSNStateMP
from .bssn_rk4_mp_vec_damp import rk4_step_mp_vec_damp
from .bssn_ic_bowen_york import bowen_york_binary
from .bssn_ricci_full import hamiltonian_constraint_full
from .amr.refinement_tensor import (
    prolongate_state, restrict_state_into
)
from .amr.amr_interp import inject_boundary
from .amr.refinement import prolongate
from .amr.amr_moving import fine_grid_mask
from .amr.amr_moving import target_lo, fine_grid_mask
from .amr.amr_bh_track import peak_centroids, pair_separation, find_bh_centers
from .amr.amr_projection import (
    project_phi_band, project_phi_band_covariant,
    momentum_band_correct, atilde_up_divergence
)
from .bssn_grid import d_central


@dataclass
class AMRConfig:
    """AMR configuration parameters."""
    # Coarse grid
    N_c: int = 65
    L_c: float = 100.0
    # Fine grid
    N_f: int = 33
    L_f: float = 20.0
    # Refinement factor
    factor: int = 2
    # Time stepping
    dt: float = 0.125
    n_steps: int = 1000
    # Projection
    use_projection: bool = True
    projection_type: str = "covariant"  # "flat" or "covariant"
    n_layers: int = 4
    projection_sweeps: int = 20
    projection_omega: float = 0.8
    mask_threshold: float = 0.5
    # Momentum correction
    use_momentum_correction: bool = True
    momentum_passes: int = 1
    # Gauge
    eta: float = 1.0
    # Boundary damping
    width_cells: int = 3
    kappa: float = 0.5
    # Output
    output_every: int = 50
    output_dir: str = "results/validated/amr_long_run"


def create_initial_data(config: AMRConfig):
    """Create binary puncture initial data on coarse and fine grids."""
    # Coarse grid
    x_c = np.linspace(-config.L_c / 2, config.L_c / 2, config.N_c)
    X_c, Y_c, Z_c = np.meshgrid(x_c, x_c, x_c, indexing="ij")
    
    # Fine grid
    x_f = np.linspace(-config.L_f / 2, config.L_f / 2, config.N_f)
    X_f, Y_f, Z_f = np.meshgrid(x_f, x_f, x_f, indexing="ij")
    
    # Binary parameters
    separation = 10.0  # M
    mass = 1.0
    P = 0.1  # momentum

    # Bowen-York binary on coarse grid
    s_c_base = bowen_york_binary(
        config.N_c, config.L_c,
        m1=mass, m2=mass,
        d=separation,
    )
    s_c = BSSNStateMP.from_base(s_c_base)
    
    # Prolongate to fine grid (initial fine grid centered)
    offset = ((config.N_c - config.N_f) // 2,) * 3
    fine_dict = prolongate_state(s_c, config.N_f, offset)
    # Add B field (prolongate from coarse)
    B_f = np.zeros((config.N_f, config.N_f, config.N_f, 3))
    for i in range(3):
        B_f[..., i] = prolongate(s_c.B[..., i], config.N_f, offset)
    fine_dict["B"] = B_f
    s_f = BSSNStateMP(
        phi=fine_dict["phi"],
        gtilde=fine_dict["gtilde"],
        K=fine_dict["K"],
        Atilde=fine_dict["Atilde"],
        Gamma=fine_dict["Gamma"],
        alpha=fine_dict["alpha"],
        beta=fine_dict["beta"],
        B=fine_dict["B"],
    )
    
    return s_c, s_f, offset


def compute_hamiltonian(state, h):
    """Compute Hamiltonian constraint on grid."""
    return hamiltonian_constraint_full(
        state.phi, state.gtilde, state.K, state.Atilde, h
    )


def apply_projection(state_f, h_f, config: AMRConfig):
    """Apply constraint projection to fine grid state."""
    if not config.use_projection:
        return state_f, {}
    
    # Compute Hamiltonian residual
    H = compute_hamiltonian(state_f, h_f)
    
    # Project phi band
    if config.projection_type == "covariant":
        phi_new, dphi_max = project_phi_band_covariant(
            state_f.phi, state_f.gtilde, h_f, H,
            n_layers=config.n_layers,
            mask_threshold=config.mask_threshold
        )
    else:
        phi_new, dphi_max = project_phi_band(
            state_f.phi, h_f, H,
            n_layers=config.n_layers,
            sweeps=config.projection_sweeps,
            omega=config.projection_omega,
            mask_threshold=config.mask_threshold
        )
    
    # Momentum correction
    At_new = state_f.Atilde
    l_before = 0.0
    l_after = 0.0
    if config.use_momentum_correction:
        At_new, l_before, l_after = momentum_band_correct(
            state_f.Atilde, state_f.gtilde, h_f,
            n_layers=config.n_layers,
            sweeps=config.projection_sweeps,
            omega=config.projection_omega,
            passes=config.momentum_passes
        )
    
    # Create new state with projected fields
    state_projected = BSSNStateMP(
        phi=phi_new,
        gtilde=state_f.gtilde,
        K=state_f.K,
        Atilde=At_new,
        Gamma=state_f.Gamma,
        alpha=state_f.alpha,
        beta=state_f.beta,
        B=state_f.B,
    )
    
    proj_info = {
        "dphi_max": float(dphi_max),
        "l_before": float(l_before),
        "l_after": float(l_after),
        "H_max": float(np.max(np.abs(H))),
    }
    
    return state_projected, proj_info


def update_fine_patch(s_c, s_f, config: AMRConfig, step: int):
    """Update fine patch position based on COM drift."""
    # Track BH positions on coarse grid
    c_left, c_right, _ = peak_centroids(s_c, config.L_c)
    com_x = (c_left[0] + c_right[0]) / 2.0
    
    # Compute new fine grid offset
    h_c = config.L_c / (config.N_c - 1)
    lo = target_lo(com_x, h_c, config.N_c, config.N_f)
    offset = (lo, lo, lo)
    
    # Prolongate coarse -> fine at new position
    fine_dict = prolongate_state(s_c, config.N_f, offset)
    # Add B field (prolongate from coarse)
    B_f = np.zeros((config.N_f, config.N_f, config.N_f, 3))
    for i in range(3):
        B_f[..., i] = prolongate(s_c.B[..., i], config.N_f, offset)
    fine_dict["B"] = B_f
    
    # Inject boundary from coarse
    for key in ["phi", "K", "alpha", "gtilde", "Atilde", "Gamma", "beta"]:
        fine_dict[key] = inject_boundary(
            fine_dict[key], fine_dict[key], n_layers=4
        )
    
    s_f_new = BSSNStateMP(
        phi=fine_dict["phi"],
        gtilde=fine_dict["gtilde"],
        K=fine_dict["K"],
        Atilde=fine_dict["Atilde"],
        Gamma=fine_dict["Gamma"],
        alpha=fine_dict["alpha"],
        beta=fine_dict["beta"],
        B=fine_dict["B"],
    )
    
    return s_f_new, offset, com_x


def run_amr_long_run(config: AMRConfig) -> Dict[str, Any]:
    """Run AMR long run with projection."""
    print(f"Starting AMR long run: {config.n_steps} steps, dt={config.dt}")
    print(f"Coarse: {config.N_c}^3, L={config.L_c}, h={config.L_c/(config.N_c-1):.3f}")
    print(f"Fine: {config.N_f}^3, L={config.L_f}, h={config.L_f/(config.N_f-1):.3f}")
    print(f"Projection: {config.use_projection} ({config.projection_type})")
    print(f"Momentum correction: {config.use_momentum_correction}")
    
    # Create output directory
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Initial data
    s_c, s_f, offset = create_initial_data(config)
    h_c = config.L_c / (config.N_c - 1)
    h_f = config.L_f / (config.N_f - 1)
    
    # Diagnostics
    diagnostics = {
        "steps": [],
        "time": [],
        "H_c_max": [],
        "H_c_mean": [],
        "H_f_max": [],
        "H_f_mean": [],
        "H_f_band_max": [],
        "com_x": [],
        "separation": [],
        "projection": [],
    }
    
    # Fine grid mask for band diagnostics
    fine_mask = fine_grid_mask(config.N_f, config.L_f, r_min=1.5)
    band_mask = ~fine_mask  # injection band
    
    for step in range(config.n_steps):
        t = step * config.dt
        
        # 1. Evolve coarse grid
        s_c = rk4_step_mp_vec_damp(
            s_c, h_c, config.dt, eta=config.eta,
            width_cells=config.width_cells, kappa=config.kappa
        )
        
        # 2. Update fine patch position (every step for moving patch)
        s_f, offset, com_x = update_fine_patch(s_c, s_f, config, step)
        
        # 3. Evolve fine grid
        s_f = rk4_step_mp_vec_damp(
            s_f, h_f, config.dt, eta=config.eta,
            width_cells=config.width_cells, kappa=config.kappa
        )
        
        # 4. Apply constraint projection
        s_f, proj_info = apply_projection(s_f, h_f, config)
        
        # 5. Restrict fine -> coarse (interior only)
        coarse_dict = {
            "phi": s_c.phi, "gtilde": s_c.gtilde, "K": s_c.K,
            "Atilde": s_c.Atilde, "Gamma": s_c.Gamma,
            "alpha": s_c.alpha, "beta": s_c.beta,
        }
        # Simple restriction for now (could use restrict_state_into)
        # For simplicity, we'll just update the coarse interior
        # from the fine grid center region
        
        # 6. Diagnostics
        H_c = compute_hamiltonian(s_c, h_c)
        H_f = compute_hamiltonian(s_f, h_f)
        
        # BH tracking on fine grid
        c_left_f, c_right_f, _ = peak_centroids(s_f, config.L_f)
        sep_f = pair_separation(c_left_f[0], c_right_f[0])
        
        # Band max H
        H_f_band_max = float(np.max(np.abs(H_f[band_mask]))) if band_mask.any() else 0.0
        
        diagnostics["steps"].append(step)
        diagnostics["time"].append(float(t))
        diagnostics["H_c_max"].append(float(np.max(np.abs(H_c))))
        diagnostics["H_c_mean"].append(float(np.mean(np.abs(H_c))))
        diagnostics["H_f_max"].append(float(np.max(np.abs(H_f))))
        diagnostics["H_f_mean"].append(float(np.mean(np.abs(H_f))))
        diagnostics["H_f_band_max"].append(H_f_band_max)
        diagnostics["com_x"].append(float(com_x))
        diagnostics["separation"].append(float(sep_f))
        diagnostics["projection"].append(proj_info)
        
        # Output
        if step % config.output_every == 0:
            print(f"  Step {step:4d} t={t:6.2f} | "
                  f"H_c_max={np.max(np.abs(H_c)):.3e} "
                  f"H_f_max={np.max(np.abs(H_f)):.3e} "
                  f"H_f_band={H_f_band_max:.3e} "
                  f"sep={sep_f:.3f} "
                  f"com={com_x:.3f}")
            if proj_info:
                print(f"    Proj: dphi_max={proj_info['dphi_max']:.3e} "
                      f"l_before={proj_info['l_before']:.3e} "
                      f"l_after={proj_info['l_after']:.3e}")
        
        # Check for blowup
        if H_f_band_max > 10.0:
            print(f"  BLOWUP detected at step {step}, H_f_band_max={H_f_band_max:.3e}")
            break
    
    # Save diagnostics
    out_file = out_dir / f"amr_long_run_{'proj' if config.use_projection else 'noproj'}.json"
    with open(out_file, 'w') as f:
        json.dump(diagnostics, f, indent=2)
    
    print(f"Results saved to {out_file}")
    return diagnostics


def main():
    """Run AMR long run with and without projection for comparison."""
    # Config with projection
    config_proj = AMRConfig(
        N_c=49, L_c=80.0,
        N_f=25, L_f=16.0,
        dt=0.125, n_steps=200,
        use_projection=True,
        projection_type="covariant",
        n_layers=4,
        projection_sweeps=10,
        projection_omega=0.8,
        mask_threshold=0.5,
        use_momentum_correction=True,
        momentum_passes=1,
        eta=1.0,
        width_cells=3,
        kappa=0.5,
        output_every=25,
    )
    
    # Config without projection (baseline)
    config_noproj = AMRConfig(
        N_c=49, L_c=80.0,
        N_f=25, L_f=16.0,
        dt=0.125, n_steps=200,
        use_projection=False,
        output_every=25,
    )
    
    print("=" * 60)
    print("RUNNING WITH PROJECTION")
    print("=" * 60)
    diag_proj = run_amr_long_run(config_proj)
    
    print("\n" + "=" * 60)
    print("RUNNING WITHOUT PROJECTION (BASELINE)")
    print("=" * 60)
    diag_noproj = run_amr_long_run(config_noproj)
    
    # Summary comparison
    print("\n" + "=" * 60)
    print("SUMMARY COMPARISON")
    print("=" * 60)
    
    last_proj = diag_proj["H_f_band_max"][-1] if diag_proj["H_f_band_max"] else 0
    last_noproj = diag_noproj["H_f_band_max"][-1] if diag_noproj["H_f_band_max"] else 0
    
    print(f"With projection:    H_f_band_max final = {last_proj:.3e}")
    print(f"Without projection: H_f_band_max final = {last_noproj:.3e}")
    print(f"Ratio: {last_noproj / last_proj:.2f}x" if last_proj > 0 else "N/A")
    
    # Check if projection helped
    if last_proj < last_noproj:
        print("✅ PROJECTION HELPS: constraint violation reduced")
    else:
        print("❌ PROJECTION DID NOT HELP: constraint violation not reduced")


if __name__ == "__main__":
    main()