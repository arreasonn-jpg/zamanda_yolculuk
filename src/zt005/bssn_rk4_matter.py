"""RK4 driver coupling scalar field to BSSN metric (v4.0-alpha1)."""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_rhs_matter import bssn_rhs_matter
from .bssn_grid import d_central


def _matter_rhs(phi_field, Pi_field, alpha, beta, h, V_func=None):
    """Scalar field evolution:
      d_t phi = alpha Pi + beta^i d_i phi
      d_t Pi  = alpha K Pi + alpha Lap(phi) + beta^i d_i Pi
                + d_i alpha d^i phi - alpha dV/dphi
    where Lap is flat laplacian (approx; conformal factor corrections
    omitted at this level).
    """
    dphi = np.stack(d_central(phi_field, h), axis=-1)
    dPi  = np.stack(d_central(Pi_field, h), axis=-1)
    dalpha = np.stack(d_central(alpha, h), axis=-1)

    lap = np.zeros_like(phi_field)
    for ax in range(3):
        dd = d_central(dphi[..., ax], h)
        lap += dd[ax]

    if V_func is None:
        dV = 0.0
    else:
        eps = 1e-4
        dV = np.array([(V_func(p+eps) - V_func(p-eps))/(2*eps)
                       for p in phi_field.ravel()]).reshape(phi_field.shape)

    # For now, K = 0 (K is from metric state; pass in later)
    K = np.zeros_like(alpha)
    bdot_dphi = np.einsum('...k,...k->...', beta, dphi)
    bdot_dPi  = np.einsum('...k,...k->...', beta, dPi)
    grad_alpha_grad_phi = np.einsum('...k,...k->...', dalpha, dphi)

    dphi_t = alpha * Pi_field + bdot_dphi
    dPi_t  = alpha * K * Pi_field + alpha * lap + bdot_dPi \
             + grad_alpha_grad_phi - alpha * dV
    return dphi_t, dPi_t


def rk4_step_matter(state: BSSNState, phi_field, Pi_field, h, dt,
                    V_func=None):
    """One RK4 step on (metric + scalar). Returns (new_state, phi, Pi)."""
    def metric_rhs(s):
        r, _ = bssn_rhs_matter(s.phi, s.gtilde, s.K, s.Atilde, s.alpha,
                                s.beta, phi_field, Pi_field, h, V_func)
        return r

    def matter_rhs(s, pf, pif):
        return _matter_rhs(pf, pif, s.alpha, s.beta, h, V_func)

    def add_metric(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy(),
        )

    # k1
    k1 = metric_rhs(state)
    m1 = matter_rhs(state, phi_field, Pi_field)
    s2 = add_metric(state, k1, 0.5*dt)
    phi2 = phi_field + 0.5*dt*m1[0]
    Pi2  = Pi_field  + 0.5*dt*m1[1]

    # k2
    k2 = metric_rhs(s2)
    m2 = matter_rhs(s2, phi2, Pi2)
    s3 = add_metric(state, k2, 0.5*dt)
    phi3 = phi_field + 0.5*dt*m2[0]
    Pi3  = Pi_field  + 0.5*dt*m2[1]

    # k3
    k3 = metric_rhs(s3)
    m3 = matter_rhs(s3, phi3, Pi3)
    s4 = add_metric(state, k3, dt)
    phi4 = phi_field + dt*m3[0]
    Pi4  = Pi_field  + dt*m3[1]

    # k4
    k4 = metric_rhs(s4)
    m4 = matter_rhs(s4, phi4, Pi4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy(),
    )
    new_phi = phi_field + (dt/6)*(m1[0]+2*m2[0]+2*m3[0]+m4[0])
    new_Pi  = Pi_field  + (dt/6)*(m1[1]+2*m2[1]+2*m3[1]+m4[1])
    return new_state, new_phi, new_Pi


__all__ = ["rk4_step_matter"]
