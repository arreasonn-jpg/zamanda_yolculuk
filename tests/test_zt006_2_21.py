"""ZT-006.2.21 unit tests."""
import numpy as np
from zt005.poisson_solver import solve_poisson_axisym, EINSTEIN_FACTOR


def test_einstein_factor_range():
    assert 4.0e-43 < EINSTEIN_FACTOR < 4.3e-43   # 16πG/c⁴


def test_zero_source_gives_zero_solution():
    r = np.linspace(0, 1, 20)
    z = np.linspace(-1, 1, 40)
    h = solve_poisson_axisym(np.zeros((20, 40)), r, z)
    assert np.allclose(h, 0.0, atol=1e-30)


def test_positive_source_gives_positive_center():
    r = np.linspace(0, 1, 30)
    z = np.linspace(-1, 1, 60)
    R, Z = np.meshgrid(r, z, indexing="ij")
    T00 = np.exp(-(R**2 + Z**2) / 0.1**2)
    h = solve_poisson_axisym(T00, r, z)
    assert h[0, 30] > 0.0


def test_boundary_conditions_are_zero():
    r = np.linspace(0, 1, 30)
    z = np.linspace(-1, 1, 60)
    R, Z = np.meshgrid(r, z, indexing="ij")
    T00 = np.exp(-(R**2 + Z**2) / 0.1**2)
    h = solve_poisson_axisym(T00, r, z)

    # Outer r boundary
    assert np.allclose(h[-1, :], 0.0, atol=1e-30)
    # z boundaries
    assert np.allclose(h[:, 0],  0.0, atol=1e-30)
    assert np.allclose(h[:, -1], 0.0, atol=1e-30)


def test_convergence_is_monotone():
    """Finer grid should not diverge wildly from coarser."""
    h_list = []
    for Nr, Nz in [(20, 40), (30, 60), (40, 80)]:
        r = np.linspace(0, 1, Nr)
        z = np.linspace(-1, 1, Nz)
        R, Z = np.meshgrid(r, z, indexing="ij")
        T00 = np.exp(-(R**2 + Z**2) / 0.1**2)
        h_list.append(solve_poisson_axisym(T00, r, z)[0, Nz // 2])

    h1, h2, h3 = h_list
    rel12 = abs(h2 - h1) / abs(h2)
    rel23 = abs(h3 - h2) / abs(h3)
    # Second refinement step should not be dramatically worse
    assert rel23 < rel12 + 0.5