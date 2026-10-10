"""RK4 driver with frozen gauge (alpha, beta from IC, not evolved)."""
from __future__ import annotations
from .bssn_state import BSSNState
from .bssn_rk4_beta2 import rk4_step_beta2
from .bssn_bc import sommerfeld_all


def rk4_step_frozen(state: BSSNState, h, dt):
    """One RK4 step with alpha, beta held fixed."""
    alpha0 = state.alpha.copy()
    beta0 = state.beta.copy()
    old = BSSNState(
        phi=state.phi.copy(), gtilde=state.gtilde.copy(),
        K=state.K.copy(), Atilde=state.Atilde.copy(),
        Gamma=state.Gamma.copy(), alpha=alpha0.copy(),
        beta=beta0.copy(),
    )
    new = rk4_step_beta2(state, h, dt)
    # restore gauge
    new.alpha = alpha0
    new.beta = beta0
    # BC on the evolved fields only
    d = sommerfeld_all(new, old, h, dt, c=1.0)
    return BSSNState(
        phi=d["phi"], gtilde=d["gtilde"], K=d["K"],
        Atilde=d["Atilde"], Gamma=d["Gamma"],
        alpha=alpha0, beta=beta0,
    )


__all__ = ["rk4_step_frozen"]
