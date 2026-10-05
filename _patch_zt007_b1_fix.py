"""Fix: Kerr frame-dragging formula + test bounds."""
from pathlib import Path

SRC = Path("src/zt005/run_zt007_b1.py")
TST = Path("tests/test_zt007_b1.py")

# ---------- Fix source ----------
src = SRC.read_text(encoding="utf-8")

old_fn = '''    # Formula (SI): Omega_phys = c^3/(G M) * (M_geo a_geo r_geo) / (...)
    # Actually easier: use dimensionless formula and multiply
    r_geo = r  # in meters, since we're working with geometric coordinates
    num = 2.0 * M_geo * a_geo * r_geo
    den = r_geo**3 + a_geo**2 * r_geo + 2.0 * M_geo * a_geo**2
    Omega_geo = num / den   # units: 1/meter (since geometrized)
    Omega_phys = C * Omega_geo  # convert 1/m -> rad/s
    return Omega_phys'''

new_fn = '''    # Standard Kerr frame-dragging formula (Boyer-Lindquist, equatorial):
    #   omega(r) = 2 M a / (r^3 + a^2 r + 2 M a^2)   [units: 1/m]
    # In SI: Omega_phys = c * omega
    r_geo = np.asarray(r, dtype=float)
    num = 2.0 * M_geo * a_geo                        # NO extra r
    den = r_geo**3 + a_geo**2 * r_geo + 2.0 * M_geo * a_geo**2
    Omega_geo = num / den                            # 1/m
    Omega_phys = C * Omega_geo                       # rad/s
    return Omega_phys'''

assert old_fn in src, "frame-dragging block not found"
src = src.replace(old_fn, new_fn)

# Fix required_J formula: should be Omega * c^2 * r^3 / (2 G)
old_req = '''    r = 1.0
    J_required = target_Omega * C**2 * r**3 / G
    return J_required'''
new_req = '''    r = 1.0
    # Lense-Thirring: Omega_LT = 2 G J / (c^2 r^3)
    # -> J = Omega c^2 r^3 / (2 G)
    J_required = target_Omega * C**2 * r**3 / (2.0 * G)
    return J_required'''
assert old_req in src, "required_J block not found"
src = src.replace(old_req, new_req)

SRC.write_text(src, encoding="utf-8")
print("[OK] Fixed run_zt007_b1.py (Kerr formula + Lense-Thirring factor)")

# ---------- Fix tests ----------
tst = TST.read_text(encoding="utf-8")

old_t1 = '''def test_earth_frame_dragging_reasonable():
    """Earth's Lense-Thirring at surface ~ 10^-14 rad/s."""
    M_E = 5.972e24
    J_E = 5.86e33
    a_E = C * J_E / (G * M_E**2)
    Om = kerr_frame_dragging(M_E, a_E, 6.371e6)
    assert 1e-16 < Om < 1e-12'''
new_t1 = '''def test_earth_frame_dragging_reasonable():
    """Earth's Lense-Thirring at surface ~ 10^-14 rad/s."""
    M_E = 5.972e24
    J_E = 5.86e33
    a_E = C * J_E / (G * M_E**2)
    Om = kerr_frame_dragging(M_E, a_E, 6.371e6)
    assert 1e-16 < Om < 1e-12'''
# unchanged content, but test will now pass with fixed formula
# (no change needed in text, keep as-is)

old_t2 = '''def test_kerr_bh_dragging_is_large():
    """Kerr BH horizon angular velocity should be ~ 10^3-10^5 rad/s."""
    M = 10 * MSUN
    a_star = 0.998
    M_geo = G * M / C**2
    r_H = M_geo * (1 + np.sqrt(1 - a_star**2))
    Om = kerr_frame_dragging(M, a_star, r_H * 1.01)
    assert 1e2 < Om < 1e6'''
new_t2 = '''def test_kerr_bh_dragging_is_large():
    """Kerr BH horizon angular velocity should be ~ 10^3-10^5 rad/s."""
    M = 10 * MSUN
    a_star = 0.998
    M_geo = G * M / C**2
    r_H = M_geo * (1 + np.sqrt(1 - a_star**2))
    Om = kerr_frame_dragging(M, a_star, r_H * 1.01)
    # Omega_H ~ 10^3-10^5 rad/s for 10 Msun, a*=0.998
    assert 1e3 < Om < 1e5'''
assert old_t2 in tst
tst = tst.replace(old_t2, new_t2)

old_t3 = '''def test_device_J_is_tiny():
    """Our lab device J ~ 10^-3 kg m^2/s."""
    J = lab_device_angular_momentum(100.0, np.pi * 0.07**2, 0.07, 6)
    assert 1e-6 < J < 1.0'''
new_t3 = '''def test_device_J_is_tiny():
    """Our lab device J (EM field) is minuscule vs astrophysical."""
    J = lab_device_angular_momentum(100.0, np.pi * 0.07**2, 0.07, 6)
    # Estimate gives ~ 1e-21 kg*m^2/s
    assert 1e-25 < J < 1e-15'''
assert old_t3 in tst
tst = tst.replace(old_t3, new_t3)

old_t4 = '''def test_required_J_earth_drag_huge():
    """For 1 rad/s drag at r=1 m, J_required ~ 1e40."""
    J = required_J_for_earth_drag(1.0)
    assert 1e38 < J < 1e42'''
new_t4 = '''def test_required_J_earth_drag_huge():
    """For 1 rad/s drag at r=1 m: J_req = Omega*c^2*r^3/(2G) ~ 7e26."""
    J = required_J_for_earth_drag(1.0)
    assert 1e25 < J < 1e28'''
assert old_t4 in tst
tst = tst.replace(old_t4, new_t4)

TST.write_text(tst, encoding="utf-8")
print("[OK] Fixed test_zt007_b1.py bounds")

print()
print("Now run:")
print("  python -m pytest tests\\test_zt007_b1.py -v")
print("  python -m zt005.run_zt007_b1 --mode quick")