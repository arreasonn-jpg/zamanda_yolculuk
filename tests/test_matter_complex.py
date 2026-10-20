import numpy as np
import pytest
from zt005.matter.matter_complex import complex_sources, charge_density
from zt005.bssn_state import minkowski_state


def test_vacuum_zero():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    src = complex_sources(Z, Z, Z, Z, s.gtilde, s.phi, h=0.1)
    assert np.allclose(src["rho"], 0.0)
    assert np.allclose(src["S_i"], 0.0)


def test_real_reduces_to_scalar():
    """If Psi_I = Pi_I = 0, complex should match real scalar formula."""
    s = minkowski_state(5)
    Psi_R = np.ones((5,5,5)) * 0.3
    Pi_R = np.ones((5,5,5)) * 0.1
    Z = np.zeros((5,5,5))
    V_func = lambda p2: 0.5 * p2  # V = 1/2 |Psi|^2
    src = complex_sources(Psi_R, Z, Pi_R, Z, s.gtilde, s.phi, h=0.1,
                          V_func=V_func)
    # Real scalar: rho = Pi^2 + grad^2 + V
    # Here: Pi_R = 0.1, grad = 0, V = 1/2 * 0.09 = 0.045
    # rho = 0.01 + 0 + 0.045 = 0.055
    assert np.allclose(src["rho"], 0.055, atol=1e-10)


def test_charge_zero_for_real():
    Psi_R = np.ones((5,5,5)) * 0.3
    Z = np.zeros((5,5,5))
    rho_Q = charge_density(Psi_R, Z, Z, Z)
    assert np.allclose(rho_Q, 0.0)


def test_charge_nonzero_for_complex():
    Psi_R = np.ones((5,5,5)) * 0.3
    Psi_I = np.ones((5,5,5)) * 0.2
    Pi_R = np.zeros((5,5,5))
    Pi_I = np.ones((5,5,5)) * 0.5
    rho_Q = charge_density(Psi_R, Psi_I, Pi_R, Pi_I)
    # 2(Psi_R Pi_I - Psi_I Pi_R) = 2(0.3*0.5 - 0.2*0) = 0.3
    assert np.allclose(rho_Q, 0.3)


def test_phase_rotation_invariance():
    """Physical rho is invariant under (Psi_R, Psi_I) rotation."""
    s = minkowski_state(5)
    Psi_R = np.ones((5,5,5)) * 0.3
    Psi_I = np.ones((5,5,5)) * 0.4
    Pi_R = np.zeros((5,5,5))
    Pi_I = np.zeros((5,5,5))
    src1 = complex_sources(Psi_R, Psi_I, Pi_R, Pi_I, s.gtilde, s.phi, h=0.1)
    # Rotate by pi/4
    c, sn = np.cos(np.pi/4), np.sin(np.pi/4)
    Psi_R2 = c*Psi_R - sn*Psi_I
    Psi_I2 = sn*Psi_R + c*Psi_I
    src2 = complex_sources(Psi_R2, Psi_I2, Pi_R, Pi_I, s.gtilde, s.phi, h=0.1)
    assert np.allclose(src1["rho"], src2["rho"])


def test_shapes():
    s = minkowski_state(7)
    Z = np.zeros((7,7,7))
    src = complex_sources(Z, Z, Z, Z, s.gtilde, s.phi, h=0.1)
    assert src["rho"].shape == (7,7,7)
    assert src["S_i"].shape == (7,7,7,3)
    assert src["S_ij"].shape == (7,7,7,3,3)
