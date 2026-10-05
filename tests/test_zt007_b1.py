
"""ZT-007B.1 tests."""
import numpy as np
from zt005.run_zt007_b1 import (
    kerr_frame_dragging, compare_sources,
    lab_device_angular_momentum, required_J_for_earth_drag,
)
from zt005.run_zt007_b1 import G, C, MSUN


def test_earth_frame_dragging_reasonable():
    """Earth's Lense-Thirring at surface ~ 10^-14 rad/s."""
    M_E = 5.972e24
    J_E = 5.86e33
    a_E = C * J_E / (G * M_E**2)
    Om = kerr_frame_dragging(M_E, a_E, 6.371e6)
    assert 1e-16 < Om < 1e-12


def test_kerr_bh_dragging_is_large():
    """Kerr BH horizon angular velocity should be ~ 10^3-10^5 rad/s."""
    M = 10 * MSUN
    a_star = 0.998
    M_geo = G * M / C**2
    r_H = M_geo * (1 + np.sqrt(1 - a_star**2))
    Om = kerr_frame_dragging(M, a_star, r_H * 1.01)
    # Omega_H ~ 10^3-10^5 rad/s for 10 Msun, a*=0.998
    assert 1e3 < Om < 1e5


def test_sources_comparison_complete():
    s = compare_sources()
    assert "Earth" in s
    assert "NeutronStar" in s
    assert "KerrBH_10Msun" in s


def test_device_J_is_tiny():
    """Our lab device J (EM field) is minuscule vs astrophysical."""
    J = lab_device_angular_momentum(100.0, np.pi * 0.07**2, 0.07, 6)
    # Estimate gives ~ 1e-21 kg*m^2/s
    assert 1e-25 < J < 1e-15


def test_required_J_earth_drag_huge():
    """For 1 rad/s drag at r=1 m: J_req = Omega*c^2*r^3/(2G) ~ 7e26."""
    J = required_J_for_earth_drag(1.0)
    assert 1e25 < J < 1e28


def test_device_to_target_ratio():
    """Device J is >20 orders below Earth J."""
    J_dev = lab_device_angular_momentum(100.0, np.pi*0.07**2, 0.07, 6)
    J_earth = 5.86e33
    assert J_dev / J_earth < 1e-30
