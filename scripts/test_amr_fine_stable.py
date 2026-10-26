#!/usr/bin/env python3
"""
AMR test with larger fine grid domain and higher resolution for stability.
"""
import numpy as np
from zt005.amr_long_run import run_amr_long_run, AMRConfig

# Test with larger fine grid domain and higher resolution
config = AMRConfig(
    # Coarse grid - stable configuration
    N_c=17, L_c=160.0,
    # Fine grid - larger domain, higher resolution
    N_f=33, L_f=120.0,
    # Time stepping
    dt=0.0625, n_steps=50,
    # Damping
    kappa=2.0,
    # Projection
    use_projection=True,
    projection_type="covariant",
    n_layers=4,
    projection_sweeps=20,
    projection_omega=0.8,
    mask_threshold=0.5,
    use_momentum_correction=True,
    momentum_passes=1,
    # Gauge
    eta=1.0,
    # Boundary damping
    width_cells=3,
    # Output
    output_every=5,
)

print("=" * 60)
print("AMR Test: Larger fine grid (N_f=33, L_f=120)")
print("=" * 60)
print(f"Coarse: N_c={config.N_c}, L_c={config.L_c}, h_c={config.L_c/config.N_c:.3f}")
print(f"Fine:   N_f={config.N_f}, L_f={config.L_f}, h_f={config.L_f/config.N_f:.3f}")
print(f"dt={config.dt}, steps={config.n_steps}, kappa={config.kappa}")
print(f"Projection: phi={config.use_projection}, momentum={config.use_momentum_correction}")
print(f"  operator={config.projection_type}, passes={config.projection_sweeps}")
print(f"  n_layers={config.n_layers}, mask_threshold={config.mask_threshold}")
print("=" * 60)

run_amr_long_run(config)