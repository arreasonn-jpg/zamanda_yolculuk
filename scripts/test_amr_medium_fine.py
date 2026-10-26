#!/usr/bin/env python3
"""
AMR test with medium fine grid - balance of resolution and speed.
"""
import numpy as np
from zt005.amr_long_run import run_amr_long_run, AMRConfig

# Medium fine grid test
config = AMRConfig(
    N_c=17, L_c=120.0,
    N_f=25, L_f=60.0,
    dt=0.0625, n_steps=25,
    kappa=2.0,
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
    output_every=5,
)

print("=" * 60)
print("AMR Test: Medium fine grid (N_f=25, L_f=60) - 25 steps")
print("=" * 60)
print(f"Coarse: N_c={config.N_c}, L_c={config.L_c}, h_c={config.L_c/config.N_c:.3f}")
print(f"Fine:   N_f={config.N_f}, L_f={config.L_f}, h_f={config.L_f/config.N_f:.3f}")
print(f"dt={config.dt}, steps={config.n_steps}, kappa={config.kappa}")
print("=" * 60)

run_amr_long_run(config)