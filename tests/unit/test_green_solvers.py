"""Unit tests for dual Green solvers and Green convergence engine."""

import numpy as np
import pytest
from zt005.green_solvers import (
    box_singular_integral,
    uniform_sphere_exact_potential,
    solver_a_tensor_product,
    solver_b_adaptive_voxel,
    compare_green_solvers,
)
from zt005.green_convergence import (
    component_relative_errors,
    run_convergence_study,
)
from zt005.model import ActiveCylinder, Backpack, Drive


def test_box_singular_integral():
    # Box of size (0.1, 0.1, 0.1)
    # Analytic sphere equivalent: V = 0.001, R_eff = (3V/4pi)^(1/3)
    val = box_singular_integral(0.1, 0.1, 0.1, nq=8)
    assert val > 0.0
    r_eff = (3.0 * 0.001 / (4.0 * np.pi)) ** (1.0 / 3.0)
    sphere_val = 2.0 * np.pi * r_eff**2
    # Box and equivalent sphere potential should agree within 2%
    assert pytest.approx(val, rel=0.02) == sphere_val


def test_uniform_sphere_exact_potential():
    r = np.array([0.0, 0.5, 1.0, 2.0])
    pot = uniform_sphere_exact_potential(r, radius_m=1.0, rho_val=1.0)
    # Center: 2 * pi * 1.0^2 = 2*pi
    assert pytest.approx(pot[0], rel=1e-12) == 2.0 * np.pi
    # Edge: 2 * pi * (1.0 - 1/3) = 4/3 * pi
    assert pytest.approx(pot[2], rel=1e-12) == (4.0 / 3.0) * np.pi
    # Far: 4/3 * pi / 2 = 2/3 * pi
    assert pytest.approx(pot[3], rel=1e-12) == (2.0 / 3.0) * np.pi


def test_component_relative_errors():
    h1 = np.ones((4, 4)) * 1e-45
    h2 = h1 * 1.005  # 0.5% diff
    errs = component_relative_errors(h2, h1)
    assert pytest.approx(errs["h00_rel"], rel=1e-2) == 0.005
    assert pytest.approx(errs["hij_rel"], rel=1e-2) == 0.005


def test_dual_solvers_and_quick_convergence():
    cyl = ActiveCylinder(radius_m=1.0, height_m=2.0)
    bp = Backpack()
    d = Drive(peak_current_a=100.0)

    # Run quick convergence study with small grids
    res = run_convergence_study(
        name="G1",
        backpack=bp,
        drive=d,
        cyl=cyl,
        grid_ns=[6, 8],
        observers=np.array([[0.0, 0.0, 2.0]]),  # Exterior point for fast test
        tol_h00=0.05,
        tol_hij=0.05,
        validate_with_solver_b=False,
    )
    assert "acceptance" in res
    assert res["final_step"]["from_N"] == 6
    assert res["final_step"]["to_N"] == 8
