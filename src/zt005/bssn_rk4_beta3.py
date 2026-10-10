"""RK4 driver with shift + Sommerfeld BC (v2.0-beta3)."""
from __future__ import annotations
from .bssn_state import BSSNState
from .bssn_rk4_beta2 import rk4_step_beta2
from .bssn_bc import sommerfeld_all


def _apply_bc(state_new, state_old, h, dt):
    d = sommerfeld_all(state_new, state_old, h, dt, c=1.0)
    return BSSNState(
        phi=d["phi"], gtilde=d["gtilde"], K=d["K"],
        Atilde=d["Atilde"], Gamma=d["Gamma"],
        alpha=d["alpha"], beta=d["beta"],
    )


def rk4_step_beta3(state: BSSNState, h, dt):
    old = BSSNState(
        phi=state.phi.copy(), gtilde=state.gtilde.copy(),
        K=state.K.copy(), Atilde=state.Atilde.copy(),
        Gamma=state.Gamma.copy(), alpha=state.alpha.copy(),
        beta=state.beta.copy(),
    )
    new = rk4_step_beta2(state, h, dt)
    return _apply_bc(new, old, h, dt)


__all__ = ["rk4_step_beta3"]
