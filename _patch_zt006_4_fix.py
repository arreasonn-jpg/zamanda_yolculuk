"""Mini-patch: fix GPS benchmark + test."""
from pathlib import Path

SRC = Path("src/zt005/run_zt006_4.py")
TST = Path("tests/test_zt006_4.py")

# Fix the benchmark dict
src = SRC.read_text(encoding="utf-8")
old = '''    return {
        "GPS satellites (net)":          4.5e-7 * SECONDS_PER_YEAR,
        "Earth surface (relative to infinity)": 7.0e-10 * SECONDS_PER_YEAR,
        "Sun surface (relative to infinity)":  1.0e-6 * SECONDS_PER_YEAR,
        "Neutron star surface":         1.5e-1 * SECONDS_PER_YEAR,
        "Stellar black hole (10 Msun, 10 km)":  3.0e7 * SECONDS_PER_YEAR,
        "Target (1 s / year)":          1.0,
    }'''
new = '''    return {
        # 38.6 us/day net drift
        "GPS satellites (net)":          38.6e-6 * 365.25,
        # ~22 ms/year
        "Earth surface (vs infinity)":   7.0e-10 * SECONDS_PER_YEAR,
        # Sun surface ~ 2.1e-6 fraction
        "Sun surface (vs infinity)":     2.1e-6 * SECONDS_PER_YEAR,
        # Neutron star surface ~ 0.15 fraction
        "Neutron star surface":          0.15 * SECONDS_PER_YEAR,
        # Order-of-magnitude
        "Stellar BH (10 Msun, r~3r_s)":  0.3 * SECONDS_PER_YEAR,
        "Target (1 s / year)":           1.0,
    }'''
assert old in src, "benchmark block not found"
SRC.write_text(src.replace(old, new), encoding="utf-8")
print("[OK] Fixed benchmark dict")

# Fix the test
tst = TST.read_text(encoding="utf-8")
old_t = '''def test_gps_dilation_is_nanoseconds_per_day():
    b = comparison_benchmarks()
    gps_per_day = b["GPS satellites (net)"] / 365.0
    assert 1e-10 < gps_per_day < 1e-6'''
new_t = '''def test_gps_dilation_is_microseconds_per_day():
    b = comparison_benchmarks()
    gps_per_day = b["GPS satellites (net)"] / 365.25
    # 38.6 us/day
    assert 3e-5 < gps_per_day < 5e-5'''
assert old_t in tst, "test block not found"
TST.write_text(tst.replace(old_t, new_t), encoding="utf-8")
print("[OK] Fixed test")

print()
print("Now run:")
print("  python -m pytest tests\\test_zt006_4.py -v")
print("  python -m zt005.run_zt006_4 --mode quick")