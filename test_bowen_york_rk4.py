#!/usr/bin/env python3
"""Test Bowen-York RK4 step on fine grid"""
import numpy as np
from src.zt005.bssn_ic_bowen_york import bowen_york_binary
from src.zt005.bssn_rk4_grid import rk4_step_grid
from src.zt005.bssn_amr import AMRGrid, prolongate, restrict, inject_boundary
from src.zt005.bssn_constraints import hamiltonian_constraint_grid

# Small test
N_c, L_c = 17, 24
N_f, L_f = 13, 12
dt = 0.03125
kappa = 2.0

print("Creating coarse grid...")
coarse = AMRGrid(N_c, L_c, 0)
print("Creating Bowen-York...")
state_c = bowen_york_binary(coarse, separation=10.0, momentum=0.5, spin=0.0)
print("Coarse H:", np.max(np.abs(hamiltonian_constraint_grid(state_c, coarse))))

# Create fine grid
print("Creating fine grid...")
fine = AMRGrid(N_f, L_f, 1)
print("Prolongating...")
state_f = prolongate(state_c, coarse, fine)
print("Fine H:", np.max(np.abs(hamiltonian_constraint_grid(state_f, fine))))

# One RK4 step
print("Running RK4 step...")
state_f_dict = {
    "phi": state_f.phi,
    "gt": state_f.gtilde,
    "K": state_f.K,
    "At": state_f.Atilde
}
state_f_new_dict = rk4_step_grid(state_f_dict, fine.h, dt, rho_grid=None, S_grid=None, eta=kappa)
print("Fine H after RK4:", np.max(np.abs(hamiltonian_constraint_grid(state_f_new_dict, fine))))
print("Done!")