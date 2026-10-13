"""RK4 driver: Brans-Dicke + PTUM (v8.0-alpha4).

Combines:
- BSSN metric
- Brans-Dicke scalar phi_BD (with conjugate Pi_BD)
- PTUM fuel Phi + proxy psi (with conjugates)

Tests whether BD's G_eff modification shifts the PTUM QFT-OK window.
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter_brans_dicke import brans_dicke_sources
from .matter_ptum import ptum_sources
from .bssn_grid import d_central


def _lap(f, h):
    grad = np.stack(d_central(f, h), axis=-1)
    L = np.zeros_like(f)
    for ax in range(3):
        L += d_central(grad[..., ax], h)[ax]
    return L


def _bd_rhs(phi_BD, Pi_BD, alpha, K, h):
    return alpha * Pi_BD, alpha * K * Pi_BD + alpha * _lap(phi_BD, h)


def _ptum_rhs(Phi, Pi_Phi, psi, Pi_psi, alpha, K, h,
              V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    from .ptum_field import V_fuel, m_eff_squared
    eps = 1e-5
    dV_f = (V_fuel(Phi+eps, V0=V0, lam=lam) -
            V_fuel(Phi-eps, V0=V0, lam=lam)) / (2*eps)
    meff2 = m_eff_squared(Phi, m2=m2, alpha=alpha_p, Eth=Eth)
    dmeff2 = -2.0 * m2 * alpha_p * Phi / (Eth**2)
    dVp_dPhi = 0.5 * psi**2 * dmeff2
    dVp_dpsi = meff2 * psi + lam_psi * psi**3
    return (alpha * Pi_Phi,
            alpha * K * Pi_Phi + alpha * _lap(Phi, h) - alpha * (dV_f + dVp_dPhi),
            alpha * Pi_psi,
            alpha * K * Pi_psi + alpha * _lap(psi, h) - alpha * dVp_dpsi)


def _combined_source(phi_BD, Pi_BD, Phi, Pi_Phi, psi, Pi_psi,
                      gtilde, phi_conf, h, omega, ptum_params):
    src_bd = brans_dicke_sources(phi_BD, Pi_BD, gtilde, phi_conf, h, omega=omega)
    pk = dict(ptum_params)
    if "alpha_p" in pk:
        pk["alpha"] = pk.pop("alpha_p")
    src_pt = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, gtilde, phi_conf, h, **pk)
    return {
        "rho": src_bd["rho"] + src_pt["rho"],
        "S_i": src_bd["S_i"] + src_pt["S_i"],
        "S_ij": src_bd["S_ij"] + src_pt["S_ij"],
        "S": src_bd["S"] + src_pt["S"],
    }


def rk4_step_bd_ptum(state, phi_BD, Pi_BD, Phi, Pi_Phi, psi, Pi_psi,
                       h, dt, omega, ptum_params):
    def rhs_all(s, phb, Pphb, Ph, PPh, ps, Pps):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src = _combined_source(phb, Pphb, Ph, PPh, ps, Pps,
                                s.gtilde, s.phi, h, omega, ptum_params)
        e4phi = np.exp(4.0 * s.phi)
        gamma = e4phi[..., None, None] * s.gtilde
        ginv = np.linalg.inv(gamma)
        rhs["dK"] = rhs["dK"] + 4.0*np.pi*s.alpha*(src["rho"] + src["S"])
        STF = src["S_ij"] - (1.0/3.0)*src["S"][..., None, None]*gamma
        rhs["dAt"] = rhs["dAt"] + 8.0*np.pi*s.alpha[..., None, None]*\
                     np.exp(-4.0*s.phi)[..., None, None]*STF
        rhs["dGamma"] = rhs["dGamma"] - 16.0*np.pi*s.alpha[..., None]*\
                        np.einsum('...ij,...j->...i', ginv, src["S_i"])
        m_bd = _bd_rhs(phb, Pphb, s.alpha, s.K, h)
        m_pt = _ptum_rhs(Ph, PPh, ps, Pps, s.alpha, s.K, h, **ptum_params)
        return rhs, m_bd, m_pt

    def add_m(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"], gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"], Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy())

    k1, bd1, pt1 = rhs_all(state, phi_BD, Pi_BD, Phi, Pi_Phi, psi, Pi_psi)
    s2 = add_m(state, k1, 0.5*dt)
    phb2 = phi_BD + 0.5*dt*bd1[0]; Pphb2 = Pi_BD + 0.5*dt*bd1[1]
    Ph2 = Phi + 0.5*dt*pt1[0]; PPh2 = Pi_Phi + 0.5*dt*pt1[1]
    ps2 = psi + 0.5*dt*pt1[2]; Pps2 = Pi_psi + 0.5*dt*pt1[3]

    k2, bd2, pt2 = rhs_all(s2, phb2, Pphb2, Ph2, PPh2, ps2, Pps2)
    s3 = add_m(state, k2, 0.5*dt)
    phb3 = phi_BD + 0.5*dt*bd2[0]; Pphb3 = Pi_BD + 0.5*dt*bd2[1]
    Ph3 = Phi + 0.5*dt*pt2[0]; PPh3 = Pi_Phi + 0.5*dt*pt2[1]
    ps3 = psi + 0.5*dt*pt2[2]; Pps3 = Pi_psi + 0.5*dt*pt2[3]

    k3, bd3, pt3 = rhs_all(s3, phb3, Pphb3, Ph3, PPh3, ps3, Pps3)
    s4 = add_m(state, k3, dt)
    phb4 = phi_BD + dt*bd3[0]; Pphb4 = Pi_BD + dt*bd3[1]
    Ph4 = Phi + dt*pt3[0]; PPh4 = Pi_Phi + dt*pt3[1]
    ps4 = psi + dt*pt3[2]; Pps4 = Pi_psi + dt*pt3[3]

    k4, bd4, pt4 = rhs_all(s4, phb4, Pphb4, Ph4, PPh4, ps4, Pps4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy())
    new_phb = phi_BD + (dt/6)*(bd1[0]+2*bd2[0]+2*bd3[0]+bd4[0])
    new_Pphb = Pi_BD + (dt/6)*(bd1[1]+2*bd2[1]+2*bd3[1]+bd4[1])
    new_Ph = Phi + (dt/6)*(pt1[0]+2*pt2[0]+2*pt3[0]+pt4[0])
    new_PPh = Pi_Phi + (dt/6)*(pt1[1]+2*pt2[1]+2*pt3[1]+pt4[1])
    new_ps = psi + (dt/6)*(pt1[2]+2*pt2[2]+2*pt3[2]+pt4[2])
    new_Pps = Pi_psi + (dt/6)*(pt1[3]+2*pt2[3]+2*pt3[3]+pt4[3])
    return (new_state, new_phb, new_Pphb, new_Ph, new_PPh, new_ps, new_Pps)


__all__ = ["rk4_step_bd_ptum"]
