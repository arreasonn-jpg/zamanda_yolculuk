"""Fix ZT-007B.3: test bounds + docstring accuracy."""
from pathlib import Path

TST = Path("tests/test_zt007_b3.py")
SRC = Path("src/zt005/run_zt007_b3.py")

# ---- Fix tests ----
tst = TST.read_text(encoding="utf-8")

tst = tst.replace(
    """def test_max_casimir_density_nanometer():
    # At d=1nm, rho ~ -1e5 J/m^3
    rho = max_casimir_density(1e-9)
    assert 1e4 < rho < 1e6""",
    """def test_max_casimir_density_nanometer():
    # At d=1nm, rho = -pi^2 hbar c / (720 d^4) ~ -4.33e8 J/m^3
    rho = max_casimir_density(1e-9)
    assert 1e8 < rho < 1e9"""
)

tst = tst.replace(
    """def test_casimir_does_not_scale_with_throat():
    # For very tiny wormholes (Planck scale), gap shrinks
    gap_macro = wormhole_casimir_gap(1e-3, 1e-9)["gap"]
    gap_planck = wormhole_casimir_gap(1e-35, 1e-9)["gap"]
    assert gap_planck < gap_macro""",
    """def test_casimir_does_not_scale_with_throat():
    # SMALLER wormhole needs HIGHER density (rho ~ 1/b0^2)
    # So gap_planck > gap_macro
    gap_macro = wormhole_casimir_gap(1e-3, 1e-9)["gap"]
    gap_planck = wormhole_casimir_gap(1e-35, 1e-9)["gap"]
    assert gap_planck > gap_macro * 1e50"""
)

TST.write_text(tst, encoding="utf-8")
print("[OK] Fixed test bounds in test_zt007_b3.py")

# ---- Fix docstring in source ----
src = SRC.read_text(encoding="utf-8")
src = src.replace(
    'print("1 mm solucan deligi icin gereken rho ~ 1e42 J/m^3.")',
    'print("1 mm solucan deligi icin gereken rho ~ 2.9e49 J/m^3 (naive volume).")'
)
src = src.replace(
    'print(f"  1 mm wormhole gerekli  : 1e42 J/m^3")',
    'print(f"  1 mm wormhole gerekli  : 2.9e49 J/m^3")'
)
src = src.replace(
    'print(f"Casimir etkisi, 1 nm plaka araliginda ~ -1e5 J/m^3 uretir.")',
    'print(f"Casimir etkisi, 1 nm plaka araliginda ~ -4.3e8 J/m^3 uretir.")'
)
src = src.replace(
    'print("Bu kombinasyon icin R > 10^3 m gerekir (neutron star boyutunda")',
    'print("Bu kombinasyon icin R > 10^3 m gerekir (neutron star boyutunda")'
)
SRC.write_text(src, encoding="utf-8")
print("[OK] Fixed docstring numbers in run_zt007_b3.py")

print()
print("Now run:")
print("  python -m pytest tests\\test_zt007_b3.py -v")
print("  python -m zt005.run_zt007_b3 --mode quick")