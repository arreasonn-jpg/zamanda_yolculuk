#!/usr/bin/env python3
"""
AMR test with Brill-Lindquist initial data (time-symmetric, more stable).
"""
import numpy as np
from zt005.amr_long_run import run_amr_long_run, AMRConfig
from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.bssn_state_mp import BSSNStateMP
from zt005.amr.refinement_tensor import prolongate_state
from zt005.amr.refinement import prolongate

# Test with Brill-Lindquist initial data
config = AMRConfig(
    N_c=17, L_c=120.0,
    N_f=25, L_f=60.0,
    dt=0.0625, n_steps=30,
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

# Override create_initial_data to use Brill-Lindquist
def create_bl_initial_data(config: AMRConfig):
    """Create Brill-Lindquist binary puncture initial data."""
    # Coarse grid
    x_c = np.linspace(-config.L_c / 2, config.L_c / 2, config.N_c)
    X_c, Y_c, Z_c = np.meshgrid(x_c, x_c, x_c, indexing="ij")
    
    # Fine grid
    x_f = np.linspace(-config.L_f / 2, config.L_f / 2, config.N_f)
    X_f, Y_f, Z_f = np.meshgrid(x_f, x_f, x_f, indexing="ij")
    
    # Binary parameters
    separation = 10.0  # M
    mass = 1.0

    # Brill-Lindquist binary on coarse grid (time-symmetric, K=0)
    s_c_base = brill_lindquist_state(
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

# Monkey-patch the function
import zt005.amr_long_run as amr_module
amr_module.create_initial_data = create_bl_initial_data

print("=" * 60)
print("AMR Test: Brill-Lindquist initial data - 30 steps")
print("=" * 60)
print(f"Coarse: N_c={config.N_c}, L_c={config.L_c}, h_c={config.L_c/config.N_c:.3f}")
print(f"Fine:   N_f={config.N_f}, L_f={config.L_f}, h_f={config.L_f/config.N_f:.3f}")
print(f"dt={config.dt}, steps={config.n_steps}, kappa={config.kappa}")
print("=" * 60)

run_amr_long_run(config)