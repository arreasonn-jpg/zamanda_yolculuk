"""RK4 driver: metric + complex scalar + PTUM (v6.6-alpha1).

State:
- BSSN metric (7 fields)
- Complex scalar Psi = Psi_R + i Psi_I (4 real fields: Psi_R, Psi_I, Pi_R, Pi_I)
- PTUM fields (4 real: Phi, Pi_Phi, psi, Pi_psi)

Combined T_munu = T_munu(complex scalar) + T_munu(PTUM).
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter_complex import complex_sources
from .matter_ptum import ptum_sources
from .bssn_grid import d_central


def _lap(f, h):
    grad = np.stack(d_central(f, h), axis=-1)
    L = np.zeros_like(f)
    for ax in range(3):
        dd = d_central(grad[..., ax], h)
        L += dd[ax]
    return L


def _complex_rhs(PR, PI, PiR, PiI, alpha, K, h, dV_func=None):
    """Complex scalar RHS: two KG equations."""
    lapR = _lap(PR, h); lapI = _lap(PI, h)
    if dV_func is None:
        dVR = np.zeros_like(PR); dVI = np.zeros_like(PI)
    else:
        dVR, dVI = dV_func(PR, PI)
    return (alpha * PiR,
            alpha * PiI,
            alpha * K * PiR + alpha * lapR - alpha * dVR,
            alpha * K * PiI + alpha * lapI - alpha * dVI)


def _ptum_rhs(Ph, PPh, ps, Pps, alpha, K, h,
              V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    from .ptum_field import V_fuel, m_eff_squared
    eps = 1e-5
    dV_f = (V_fuel(Ph+eps, V0=V0, lam=lam) -
            V_fuel(Ph-eps, V0=V0, lam=lam)) / (2*eps)
    meff2 = m_eff_squared(Ph, m2=m2, alpha=alpha_p, Eth=Eth)
    dmeff2 = -2.0 * m2 * alpha_p * Ph / (Eth**2)
    dVp_dPhi = 0.5 * ps**2 * dmeff2
    dVp_dpsi = meff2 * ps + lam_psi * ps**3
    return (alpha * PPh,
            alpha * K * PPh + alpha * _lap(Ph, h) - alpha * (dV_f + dVp_dPhi),
            alpha * Pps,
            alpha * K * Pps + alpha * _lap(ps, h) - alpha * dVp_dpsi)


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


def rk4_step_complex_ptum(state,
                            PR, PI, PiR, PiI,
                            Ph, PPh, ps, Pps,
                            h, dt, ptum_params, dV_func=None):
    """One RK4 step on (metric + 4 complex + 4 PTUM fields)."""
    def rhs_all(s, pr, pi, pir, pii, ph, pph, pss, pps_):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src_c = complex_sources(pr, pi, pir, pii, s.gtilde, s.phi, h)
        ptum_kwargs = dict(ptum_params)
        if "alpha_p" in ptum_kwargs:
            ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")
        src_p = ptum_sources(ph, pph, pss, pps_, s.gtilde, s.phi, h,
                              **ptum_kwargs)
        src = {
            "rho": src_c["rho"] + src_p["rho"],
            "S_i": src_c["S_i"] + src_p["S_i"],
            "S_ij": src_c["S_ij"] + src_p["S_ij"],
            "S": src_c["S"] + src_p["S"],
        }
        rhs = _apply_metric_source(rhs, src, s)
        mc = _complex_rhs(pr, pi, pir, pii, s.alpha, s.K, h, dV_func)
        mp = _ptum_rhs(ph, pph, pss, pps_, s.alpha, s.K, h, **ptum_params)
        return rhs, mc, mp

    def add_metric(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy(),
        )

    k1, c1, p1 = rhs_all(state, PR, PI, PiR, PiI, Ph, PPh, ps, Pps)
    s2 = add_metric(state, k1, 0.5*dt)
    PR2 = PR+0.5*dt*c1[0]; PI2 = PI+0.5*dt*c1[1]
    PiR2 = PiR+0.5*dt*c1[2]; PiI2 = PiI+0.5*dt*c1[3]
    Ph2 = Ph+0.5*dt*p1[0]; PPh2 = PPh+0.5*dt*p1[1]
    ps2 = ps+0.5*dt*p1[2]; Pps2 = Pps+0.5*dt*p1[3]

    k2, c2, p2 = rhs_all(s2, PR2, PI2, PiR2, PiI2, Ph2, PPh2, ps2, Pps2)
    s3 = add_metric(state, k2, 0.5*dt)
    PR3 = PR+0.5*dt*c2[0]; PI3 = PI+0.5*dt*c2[1]
    PiR3 = PiR+0.5*dt*c2[2]; PiI3 = PiI+0.5*dt*c2[3]
    Ph3 = Ph+0.5*dt*p2[0]; PPh3 = PPh+0.5*dt*p2[1]
    ps3 = ps+0.5*dt*p2[2]; Pps3 = Pps+0.5*dt*p2[3]

    k3, c3, p3 = rhs_all(s3, PR3, PI3, PiR3, PiI3, Ph3, PPh3, ps3, Pps3)
    s4 = add_metric(state, k3, dt)
    PR4 = PR+dt*c3[0]; PI4 = PI+dt*c3[1]
    PiR4 = PiR+dt*c3[2]; PiI4 = PiI+dt*c3[3]
    Ph4 = Ph+dt*p3[0]; PPh4 = PPh+dt*p3[1]
    ps4 = ps+dt*p3[2]; Pps4 = Pps+dt*p3[3]

    k4, c4, p4 = rhs_all(s4, PR4, PI4, PiR4, PiI4, Ph4, PPh4, ps4, Pps4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy(),
    )
    new_PR = PR + (dt/6)*(c1[0]+2*c2[0]+2*c3[0]+c4[0])
    new_PI = PI + (dt/6)*(c1[1]+2*c2[1]+2*c3[1]+c4[1])
    new_PiR = PiR + (dt/6)*(c1[2]+2*c2[2]+2*c3[2]+c4[2])
    new_PiI = PiI + (dt/6)*(c1[3]+2*c2[3]+2*c3[3]+c4[3])
    new_Ph = Ph + (dt/6)*(p1[0]+2*p2[0]+2*p3[0]+p4[0])
    new_PPh = PPh + (dt/6)*(p1[1]+2*p2[1]+2*p3[1]+p4[1])
    new_ps = ps + (dt/6)*(p1[2]+2*p2[2]+2*p3[2]+p4[2])
    new_Pps = Pps + (dt/6)*(p1[3]+2*p2[3]+2*p3[3]+p4[3])
    return (new_state, new_PR, new_PI, new_PiR, new_PiI,
            new_Ph, new_PPh, new_ps, new_Pps)


__all__ = ["rk4_step_complex_ptum"]
