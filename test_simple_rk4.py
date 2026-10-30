#!/usr/bin/env python3
"""Simple test of Brill-Lindquist initial data and RK4"""
import numpy as np
from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.bssn_rhs_vec import bssn_rhs_vec
from zt005.bssn_constraints import hamiltonian_constraint_grid

N = 33
L = 80.0
h = L / (N - 1)

print("Creating Brill-Lindquist...")
state = brill_lindquist_state(N, L, 0.5, 0.5, 12.0)
print("State created:", state.phi.shape)

H = hamiltonian_constraint_grid(state.phi, state.gtilde, state.K, state.Atilde, h)
print("Initial H max:", np.max(np.abs(H)))

# One RK4 step
print("Running RK4 step...")
rhs = bssn_rhs_vec(
    state.phi, state.gtilde, state.K, state.Atilde,
    state.alpha, h,
    sponge=None, kappa=2.0,
    beta=state.beta, Gamma=state.Gamma
)
print("RHS computed")

state.phi += 0.03125 * rhs["dphi"]
state.K += 0.03125 * rhs["dK"]
state.gtilde += 0.03125 * rhs["dgt"]
state.Atilde += 0.03125 * rhs["dAt"]
state.Gamma += 0.03125 * rhs["dGamma"]
print("State updated")

H = hamiltonian_constraint_grid(state.phi, state.gtilde, state.K, state.Atilde, h)
print("H after 1 step:", np.max(np.abs(H)))
print("Done!")