"""Vectorized moving puncture RK4 + boundary constraint damping
(v12.7-alpha2).

After each RK4 step, apply a sponge-layer constraint damping to K
near the boundary, and multi-layer Sommerfeld to all evolved fields.
"""
from __future__ import annotations
import numpy as np
from .bssn_state_mp import BSSNStateMP
from .bssn_rk4_mp_vec import rk4_step_mp_vec
from .bssn_bc_improved import (
    sponge_3d, damp_K, multi_layer_sommerfeld,
)
from .bssn_ricci_full import hamiltonian_constraint_full

def _apply_bc(state, state_old, h, dt, N, width_cells=3,
              kappa=0.5):
    """Apply sponge + multi-layer Sommerfeld to state fields."""
    sponge = sponge_3d(N, width_cells)
    H = hamiltonian_constraint_full(
        state.phi, state.gtilde, state.K, state.Atilde, h)
    K_new = damp_K(state.K, H, sponge, kappa=kappa)
    # Multi-layer Sommerfeld on K, phi, gtilde diag
    K_new = multi_layer_sommerfeld(K_new, state_old.K, h, dt, width_cells)
    phi_new = multi_layer_sommerfeld(state.phi, state_old.phi, h, dt, width_cells)
    return BSSNStateMP(
        phi=phi_new,
        gtilde=state.gtilde,
        K=K_new,
        Atilde=state.Atilde,
        Gamma=state.Gamma,
        alpha=state.alpha,
        beta=state.beta,
        B=state.B,
    )

def rk4_step_mp_vec_damp(state, h, dt, eta=1.0,
                          width_cells=3, kappa=0.5):
    """One RK4 step + boundary constraint damping."""
    state_old = BSSNStateMP(
        phi=state.phi.copy(), gtilde=state.gtilde.copy(),
        K=state.K.copy(), Atilde=state.Atilde.copy(),
        Gamma=state.Gamma.copy(), alpha=state.alpha.copy(),
        beta=state.beta.copy(), B=state.B.copy(),
    )
    new = rk4_step_mp_vec(state, h, dt, eta=eta)
    N = state.phi.shape[0]
    return _apply_bc(new, state_old, h, dt, N,
                     width_cells=width_cells, kappa=kappa)

__all__ = ["rk4_step_mp_vec_damp", "_apply_bc"]
