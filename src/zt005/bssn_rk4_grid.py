"""RK4 with the grid-aware full BSSN RHS."""
from __future__ import annotations
import numpy as np
from .bssn_rhs_grid import bssn_rhs_grid


def _add(s, k, dt):
    return {key: s[key] + dt * k[key] for key in ("phi", "gt", "K", "At")}


def rk4_step_grid(state, h, dt, rho_grid=None, S_grid=None, eta=0.1):
    """One RK4 step. state = dict(phi, gt, K, At)."""
    k1 = bssn_rhs_grid(state["phi"], state["gt"], state["K"], state["At"],
                       h, rho_grid, S_grid, eta)
    s2 = {"phi": state["phi"] + 0.5*dt*k1["dphi"],
          "gt":  state["gt"]  + 0.5*dt*k1["dgt"],
          "K":   state["K"]   + 0.5*dt*k1["dK"],
          "At":  state["At"]  + 0.5*dt*k1["dAt"]}
    k2 = bssn_rhs_grid(s2["phi"], s2["gt"], s2["K"], s2["At"],
                       h, rho_grid, S_grid, eta)
    s3 = {"phi": state["phi"] + 0.5*dt*k2["dphi"],
          "gt":  state["gt"]  + 0.5*dt*k2["dgt"],
          "K":   state["K"]   + 0.5*dt*k2["dK"],
          "At":  state["At"]  + 0.5*dt*k2["dAt"]}
    k3 = bssn_rhs_grid(s3["phi"], s3["gt"], s3["K"], s3["At"],
                       h, rho_grid, S_grid, eta)
    s4 = {"phi": state["phi"] + dt*k3["dphi"],
          "gt":  state["gt"]  + dt*k3["dgt"],
          "K":   state["K"]   + dt*k3["dK"],
          "At":  state["At"]  + dt*k3["dAt"]}
    k4 = bssn_rhs_grid(s4["phi"], s4["gt"], s4["K"], s4["At"],
                       h, rho_grid, S_grid, eta)
    return {
        "phi": state["phi"] + (dt/6)*(k1["dphi"] + 2*k2["dphi"] + 2*k3["dphi"] + k4["dphi"]),
        "gt":  state["gt"]  + (dt/6)*(k1["dgt"]  + 2*k2["dgt"]  + 2*k3["dgt"]  + k4["dgt"]),
        "K":   state["K"]   + (dt/6)*(k1["dK"]   + 2*k2["dK"]   + 2*k3["dK"]   + k4["dK"]),
        "At":  state["At"]  + (dt/6)*(k1["dAt"]  + 2*k2["dAt"]  + 2*k3["dAt"]  + k4["dAt"]),
    }


__all__ = ["rk4_step_grid"]
