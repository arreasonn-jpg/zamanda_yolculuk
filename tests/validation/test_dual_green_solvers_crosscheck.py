"""Validation test: cross-comparison of independent Green Solver A and Solver B."""

import numpy as np
import pytest
from zt005.manufactured_sources import GaussianBump, conserved_T
from zt005.gauss_cyl_quadrature import gauss_cylindrical_nodes
from zt005.green_solvers import (
    solver_a_tensor_product,
    solver_b_adaptive_voxel,
    compare_green_solvers,
)


def test_dual_solvers_crosscheck():
    bump = GaussianBump(A=1.0, sigma_t=1e-3, sigma_x=0.5)
    obs = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    cyl_r, cyl_h = 1.0, 2.5

    pts_a, w_a = gauss_cylindrical_nodes(cyl_r, cyl_h, 28, 56, 28)
    T_a = conserved_T(bump, np.zeros(len(pts_a)), pts_a)
    h_a = solver_a_tensor_product(obs, pts_a, w_a, T_a)

    def eval_bump(pts):
        return conserved_T(bump, np.zeros(len(pts)), pts)

    h_b = solver_b_adaptive_voxel(obs, cyl_r, cyl_h, eval_bump, nx=36, ny=36, nz=36)
    cmp = compare_green_solvers(h_a, h_b, eps_tol=0.01)

    assert cmp["pass"] is True
    assert cmp["tensor_rel_err"] < 0.01
    assert cmp["h00_rel_err"] < 0.01
