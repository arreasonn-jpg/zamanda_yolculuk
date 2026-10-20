"""RK4 driver: Kerr+SC metric + SC condensate + PTUM fields (v6.4-alpha4).

Evolves:
- BSSN metric (7 fields)
- SC condensate: phi_sc, Pi_sc (real scalar)
- PTUM: Phi, Pi_Phi, psi, Pi_psi (4 real fields)

Combined source terms: T_munu(SC) + T_munu(PTUM).
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter.matter_scalar import scalar_sources
from .matter.matter_ptum import ptum_sources
from .bssn_grid import d_central


def _lap(f, h):
    grad = np.stack(d_central(f, h), axis=-1)
    L = np.zeros_like(f)
    for ax in range(3):
        dd = d_central(grad[..., ax], h)
        L += dd[ax]
    return L


def _scalar_rhs(phi, Pi, alpha, K, h, V_func):
    """dphi = alpha Pi; dPi = alpha K Pi + alpha Lap(phi) - alpha dV/dphi."""
    eps = 1e-5
    dV = (np.array([V_func(p+eps) - V_func(p-eps)
                    for p in phi.ravel()]).reshape(phi.shape)) / (2*eps)
    return alpha * Pi, alpha * K * Pi + alpha * _lap(phi, h) - alpha * dV


def _ptum_rhs(Phi, Pi_Phi, psi, Pi_psi, alpha, K, h,
              V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    """Two-field PTUM RHS (from matter_ptum)."""
    from .matter.ptum_field import V_fuel, m_eff_squared
    eps = 1e-5
    dV_f = (V_fuel(Phi+eps, V0=V0, lam=lam) -
            V_fuel(Phi-eps, V0=V0, lam=lam)) / (2*eps)
    meff2 = m_eff_squared(Phi, m2=m2, alpha=alpha_p, Eth=Eth)
    dmeff2_dPhi = -2.0 * m2 * alpha_p * Phi / (Eth**2)
    dVp_dPhi = 0.5 * psi**2 * dmeff2_dPhi
    dVp_dpsi = meff2 * psi + lam_psi * psi**3
    return (alpha * Pi_Phi,
            alpha * K * Pi_Phi + alpha * _lap(Phi, h) - alpha * (dV_f + dVp_dPhi),
            alpha * Pi_psi,
            alpha * K * Pi_psi + alpha * _lap(psi, h) - alpha * dVp_dpsi)


def _combined_source(phi_sc, Pi_sc, Phi, Pi_Phi, psi, Pi_psi,
                      gtilde, phi_conf, h, V_sc, ptum_params):
    src_sc = scalar_sources(phi_sc, Pi_sc, gtilde, phi_conf, h, V_func=V_sc)
    ptum_kwargs = dict(ptum_params)
    if "alpha_p" in ptum_kwargs:
        ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")
    src_pt = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, gtilde, phi_conf, h,
                          **ptum_kwargs)
    return {
        "rho": src_sc["rho"] + src_pt["rho"],
        "S_i": src_sc["S_i"] + src_pt["S_i"],
        "S_ij": src_sc["S_ij"] + src_pt["S_ij"],
        "S": src_sc["S"] + src_pt["S"],
    }


def _apply_metric_source(rhs, src, state):
    e4phi = np.exp(4.0 * state.phi)
    gamma = e4phi[..., None, None] * state.gtilde
    ginv = np.linalg.inv(gamma)
    rhs["dK"] = rhs["dK"] + 4.0 * np.pi * state.alpha * (src["rho"] + src["S"])
    STF = src["S_ij"] - (1.0/3.0) * src["S"][..., None, None] * gamma
    rhs["dAt"] = rhs["dAt"] + 8.0 * np.pi * state.alpha[..., None, None] * \
                 np.exp(-4.0*state.phi)[..., None, None] * STF
    rhs["dGamma"] = rhs["dGamma"] - 16.0 * np.pi * state.alpha[..., None] * \
                    np.einsum('...ij,...j->...i', ginv, src["S_i"])
    return rhs


def rk4_step_sc_ptum(state, phi_sc, Pi_sc, Phi, Pi_Phi, psi, Pi_psi,
                      h, dt, V_sc, ptum_params):
    """One RK4 step on (metric + 6 matter fields)."""
    def rhs_all(s, psc, Psc, Ph, PPh, ps, Pps):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src = _combined_source(psc, Psc, Ph, PPh, ps, Pps,
                                s.gtilde, s.phi, h, V_sc, ptum_params)
        rhs = _apply_metric_source(rhs, src, s)
        m_sc = _scalar_rhs(psc, Psc, s.alpha, s.K, h, V_sc)
        m_pt = _ptum_rhs(Ph, PPh, ps, Pps, s.alpha, s.K, h, **ptum_params)
        return rhs, m_sc, m_pt

    def add_metric(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy(),
        )

    k1, sc1, pt1 = rhs_all(state, phi_sc, Pi_sc, Phi, Pi_Phi, psi, Pi_psi)
    s2 = add_metric(state, k1, 0.5*dt)
    psc2 = phi_sc + 0.5*dt*sc1[0]; Psc2 = Pi_sc + 0.5*dt*sc1[1]
    Ph2 = Phi + 0.5*dt*pt1[0]; PPh2 = Pi_Phi + 0.5*dt*pt1[1]
    ps2 = psi + 0.5*dt*pt1[2]; Pps2 = Pi_psi + 0.5*dt*pt1[3]

    k2, sc2, pt2 = rhs_all(s2, psc2, Psc2, Ph2, PPh2, ps2, Pps2)
    s3 = add_metric(state, k2, 0.5*dt)
    psc3 = phi_sc + 0.5*dt*sc2[0]; Psc3 = Pi_sc + 0.5*dt*sc2[1]
    Ph3 = Phi + 0.5*dt*pt2[0]; PPh3 = Pi_Phi + 0.5*dt*pt2[1]
    ps3 = psi + 0.5*dt*pt2[2]; Pps3 = Pi_psi + 0.5*dt*pt2[3]

    k3, sc3, pt3 = rhs_all(s3, psc3, Psc3, Ph3, PPh3, ps3, Pps3)
    s4 = add_metric(state, k3, dt)
    psc4 = phi_sc + dt*sc3[0]; Psc4 = Pi_sc + dt*sc3[1]
    Ph4 = Phi + dt*pt3[0]; PPh4 = Pi_Phi + dt*pt3[1]
    ps4 = psi + dt*pt3[2]; Pps4 = Pi_psi + dt*pt3[3]

    k4, sc4, pt4 = rhs_all(s4, psc4, Psc4, Ph4, PPh4, ps4, Pps4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy(),
    )
    new_phi_sc = phi_sc + (dt/6)*(sc1[0]+2*sc2[0]+2*sc3[0]+sc4[0])
    new_Pi_sc = Pi_sc + (dt/6)*(sc1[1]+2*sc2[1]+2*sc3[1]+sc4[1])
    new_Phi = Phi + (dt/6)*(pt1[0]+2*pt2[0]+2*pt3[0]+pt4[0])
    new_PiPhi = Pi_Phi + (dt/6)*(pt1[1]+2*pt2[1]+2*pt3[1]+pt4[1])
    new_psi = psi + (dt/6)*(pt1[2]+2*pt2[2]+2*pt3[2]+pt4[2])
    new_Pipsi = Pi_psi + (dt/6)*(pt1[3]+2*pt2[3]+2*pt3[3]+pt4[3])

    return (new_state, new_phi_sc, new_Pi_sc,
            new_Phi, new_PiPhi, new_psi, new_Pipsi)


__all__ = ["rk4_step_sc_ptum"]
