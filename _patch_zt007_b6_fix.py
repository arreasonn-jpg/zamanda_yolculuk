"""Fix ZT-007B.6 test bounds."""
from pathlib import Path

TST = Path("tests/test_zt007_b6.py")
c = TST.read_text(encoding="utf-8")

old1 = '''def test_lab_ctc_density_is_astronomical():
    """rho required for 1 m CTC is ~ 1e21 kg/m^3."""
    omega = godel_omega_for_ctc_at_R(1.0)
    rho = godel_density(omega)
    # ~ 1e21 kg/m^3
    assert 1e20 < rho < 1e22'''

new1 = '''def test_lab_ctc_density_is_astronomical():
    """rho required for 1 m CTC is ~ 8e25 kg/m^3."""
    omega = godel_omega_for_ctc_at_R(1.0)
    rho = godel_density(omega)
    assert 1e25 < rho < 1e27'''

old2 = '''def test_lab_ctc_mass_is_huge():
    """For 1 m, mass required ~ 1e21 kg (many Earths)."""
    omega = godel_omega_for_ctc_at_R(1.0)
    rho = godel_density(omega)
    m = mass_of_uniform_sphere(rho, 1.0)
    # ~ 4e21 kg
    assert 1e20 < m < 1e23'''

new2 = '''def test_lab_ctc_mass_is_huge():
    """For 1 m, mass required ~ 3.5e26 kg (58 Earths)."""
    omega = godel_omega_for_ctc_at_R(1.0)
    rho = godel_density(omega)
    m = mass_of_uniform_sphere(rho, 1.0)
    assert 1e25 < m < 1e28'''

assert old1 in c and old2 in c
c = c.replace(old1, new1).replace(old2, new2)
TST.write_text(c, encoding="utf-8")
print("[OK] Fixed Godel test bounds")
print()
print("Now run:")
print("  python -m pytest tests\\test_zt007_b6.py -v")