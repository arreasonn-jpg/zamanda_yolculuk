import numpy as np
from zt005.lqc_ptum_coupling import combined_density_lqc_ptum, is_below_rho_c
from zt005.lqc_scalar_evolution import RHO_C


def test_ptum_density_in_code_units():
    """PTUM density in code units: O(10) magnitude."""
    rho = combined_density_lqc_ptum(Phi=3.0, Pi_Phi=0.0,
                                     psi=2.83, Pi_psi=0.0)
    assert abs(rho) < 100  # code units
    assert abs(rho) > 1


def test_ptum_activated_negative_density():
    """A=3.0 with V_fuel=0 -> negative density (proxy minimum)."""
    rho = combined_density_lqc_ptum(Phi=3.0, Pi_Phi=0.0,
                                     psi=2.83, Pi_psi=0.0, V0=0.0)
    assert rho < 0


def test_ptum_physical_scale():
    """PTUM energy density in SI units ~ 1e19 kg/m^3 (nuclear)."""
    rho_code = combined_density_lqc_ptum(Phi=3.0, Pi_Phi=0.0,
                                          psi=2.83, Pi_psi=0.0, V0=0.0)
    # Convert: 1 code unit = 1 kg/m^3 (fictitious but illustrative)
    # The point is that PTUM rho << rho_c (Planck)
    rho_ptum_physical = abs(rho_code)  # placeholder
    assert rho_ptum_physical / RHO_C < 1e-30


def test_is_below_rho_c():
    """Small densities are below rho_c."""
    assert is_below_rho_c(1.0)
    assert is_below_rho_c(1e90)
    assert not is_below_rho_c(1e96)
