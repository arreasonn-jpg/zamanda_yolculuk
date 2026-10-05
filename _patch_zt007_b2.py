"""Auto-patch for ZT-007B.2 — Tipler cylinder CTC bound."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt007_b2.py"
TST = ROOT / "tests" / "test_zt007_b2.py"

SRC_CODE = r'''
"""ZT-007B.2 - Tipler cylinder: density-rotation bound for CTCs."""
import argparse, json
from pathlib import Path
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
C2_OVER_G = C**2 / G   # ~1.35e27 kg/m


def ctc_threshold(rho, omega, R):
    """
    Tipler CTC condition (order of magnitude, van Stockum exterior):

        r_CTC = 4 G J / c^2  must exceed R, where J = (1/2) rho omega R^4
        -> rho * omega * R^3 > c^2 / (2 G)

    Returns
    -------
    dict with:
        - 'lhs' : rho * omega * R^3
        - 'rhs' : c^2 / (2 G)
        - 'satisfied' : bool
        - 'margin' : lhs / rhs
    """
    lhs = rho * omega * R**3
    rhs = C2_OVER_G / 2.0
    return {
        "lhs": float(lhs),
        "rhs": float(rhs),
        "satisfied": bool(lhs > rhs),
        "margin": float(lhs / rhs),
    }


def critical_omega(rho, R):
    """Given rho and R, minimum omega for CTCs."""
    return C2_OVER_G / (2.0 * rho * R**3)


def critical_rho(omega, R):
    """Given omega and R, minimum density for CTCs."""
    return C2_OVER_G / (2.0 * omega * R**3)


def critical_R(rho, omega):
    """Given rho and omega, minimum cylinder radius for CTCs."""
    return (C2_OVER_G / (2.0 * rho * omega)) ** (1.0 / 3.0)


def materials_table():
    """Reference densities for context."""
    return {
        "Water":                    1.0e3,
        "Steel":                    7.9e3,
        "Lead":                     1.1e4,
        "Osmium (densest metal)":   2.26e4,
        "Earth core":               1.1e4,
        "Sun core":                 1.5e5,
        "White dwarf":              1.0e9,
        "Nuclear saturation":       2.3e17,
        "Neutron star core":        5.0e17,
        "Quark-gluon plasma":       1.0e18,
        "Planck density":           5.2e96,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 78)
    print("ZAMANDA YOLCULUK - ZT-007B.2")
    print("TIPLER CYLINDER: CTC DENSITY-ROTATION BOUND")
    print("=" * 78)
    print()
    print("Tipler condition (order of magnitude):")
    print("    rho * omega * R^3 > c^2 / (2G)")
    print(f"    c^2 / (2G) = {C2_OVER_G/2:.4e} kg/s")
    print()

    # --- Critical omega for known densities ---
    print("=" * 78)
    print("CRITICAL OMEGA (rad/s) for R = 1 m, various densities:")
    print("=" * 78)
    print(f"{'Material':<28s} {'rho (kg/m^3)':>16s} {'omega_crit (rad/s)':>22s}")
    print("-" * 78)
    R = 1.0
    for name, rho in materials_table().items():
        wc = critical_omega(rho, R)
        print(f"{name:<28s} {rho:>16.4e} {wc:>22.4e}")
    print()

    # --- Critical density for realistic omega ---
    print("=" * 78)
    print("CRITICAL DENSITY (kg/m^3) for R = 1 m, various omega:")
    print("=" * 78)
    print(f"{'omega (rad/s)':>18s}  {'Scenario':<35s} {'rho_crit (kg/m^3)':>20s}")
    print("-" * 78)
    omegas = [
        (1.0e-3, "Slow rotation"),
        (1.0,    "1 rad/s"),
        (1.0e3,  "Neutron star spin"),
        (1.0e4,  "Millisecond pulsar"),
        (1.0e6,  "Extreme lab limit"),
        (1.0e10, "Quark-gluon plasma scale"),
    ]
    for om, desc in omegas:
        rc = critical_rho(om, R)
        print(f"{om:>18.4e}  {desc:<35s} {rc:>20.4e}")
    print()

    # --- Sample scenarios ---
    print("=" * 78)
    print("SAMPLE SCENARIOS (rho, omega, R) -> CTC?")
    print("=" * 78)
    print(f"{'rho (kg/m^3)':>16s} {'omega (rad/s)':>16s} {'R (m)':>10s} "
          f"{'lhs':>14s} {'margin':>12s} {'CTC?':>8s}")
    print("-" * 78)
    scenarios = [
        (1e3,    1.0,   1.0),
        (1e4,    1.0e3, 1.0),
        (1e17,   1.0e3, 1.0),
        (1e18,   1.0e3, 1.0),
        (1e18,   1.0e6, 1.0),
        (1e18,   1.0e10, 1.0),
        (1e18,   1.0e3, 1.0e3),   # bigger cylinder
        (1e18,   1.0e6, 1.0e3),
        (1e18,   1.0e10, 1.0e3),
    ]
    for rho, om, R_s in scenarios:
        th = ctc_threshold(rho, om, R_s)
        ok = "YES" if th["satisfied"] else "no"
        print(f"{rho:>16.4e} {om:>16.4e} {R_s:>10.4e} "
              f"{th['lhs']:>14.4e} {th['margin']:>12.4e} {ok:>8s}")
    print()

    # --- Inverse: for realistic material, what R is required? ---
    print("=" * 78)
    print("INVERSE PROBLEM: Required cylinder size R_crit for real materials")
    print("(omega = 10^3 rad/s, neutron-star spin rate)")
    print("=" * 78)
    omega = 1.0e3
    for name, rho in materials_table().items():
        Rc = critical_R(rho, omega)
        flag = " (astrophysical)" if Rc < 1.0e6 else ""
        if Rc < 1.0e30:
            print(f"  {name:<28s}  R_crit = {Rc:.4e} m{flag}")
        else:
            print(f"  {name:<28s}  R_crit = {Rc:.4e} m  [absurd]")
    print()
    print("Karsilastirma: Gozlemlenebilir evren capi ~ 8.8e26 m")
    print()

    # --- Final summary ---
    print("=" * 78)
    print("SONUC")
    print("=" * 78)
    print("Tipler silindiri egzotik madde gerektirmez, fakat:")
    print()
    print("  1. Sonsuz uzunlukta olmalidir (fiziksel degil)")
    print("  2. Yuzey hizinda inanilmaz yogunluk X acisal hiz urunu ister")
    print()
    print("Nukleer yogunluk (10^18 kg/m^3) + 10^3 rad/s + R=1 m:")
    th = ctc_threshold(1e18, 1e3, 1.0)
    print(f"    margin = {th['margin']:.4e}   ->   CTC: {'YES' if th['satisfied'] else 'NO'}")
    print()
    print("Bu kombinasyon icin R > 10^3 m gerekir (neutron star boyutunda")
    print("tek bir silindir, pratikte imkansiz)")
    print()
    print("Sonuc: Tipler silindiri klasik GR'de CTC icin en 'az' egzotik")
    print("cozum olsa bile, laboratuvar olceginde UYGULANAMAZ.")
    print("=" * 78)

    out = Path("zt007_b2_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-007B.2",
        "mode": a.mode,
        "c2_over_2G": C2_OVER_G / 2.0,
        "materials": materials_table(),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: Tipler feasibility; no CTC construction.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-007B.2 tests - Tipler cylinder."""
import numpy as np
from zt005.run_zt007_b2 import (
    ctc_threshold, critical_omega, critical_rho, critical_R,
    materials_table, C2_OVER_G,
)


def test_c2_over_G_order():
    # c^2/G ~ 1.35e27 kg/m
    assert 1e27 < C2_OVER_G < 2e27


def test_ctc_threshold_satisfied_when_product_huge():
    # Huge density and rotation
    th = ctc_threshold(1e20, 1e10, 1.0)
    assert th["satisfied"] is True
    assert th["margin"] > 1.0


def test_ctc_threshold_not_satisfied_for_water():
    # Water, 1 rad/s, R=1m - hopelessly low
    th = ctc_threshold(1e3, 1.0, 1.0)
    assert th["satisfied"] is False
    assert th["margin"] < 1e-20


def test_critical_omega_for_neutron_star_density():
    # rho = 10^18 kg/m^3, R=1 m -> omega_crit = c^2/(2G rho R^3)
    wc = critical_omega(1e18, 1.0)
    expected = (2.99792458e8)**2 / (6.67430e-11) / (2 * 1e18)
    assert abs(wc - expected) / expected < 1e-6
    # Around 6.7e8 rad/s
    assert 1e8 < wc < 1e9


def test_critical_density_inverse():
    # At omega = 10^3 rad/s, R=1 m, what density?
    rc = critical_rho(1e3, 1.0)
    assert 1e23 < rc < 1e25


def test_critical_R_neutron_star():
    # Nuc density + ms spin: R_crit should be small-ish
    Rc = critical_R(1e18, 1e3)
    # Check it's the cube root of a huge number
    assert 1e2 < Rc < 1e4


def test_materials_include_neutron_star():
    m = materials_table()
    assert "Neutron star core" in m
    assert m["Neutron star core"] > 1e17


def test_increasing_R_helps_ctc():
    # Same rho, omega, larger R -> higher margin
    m1 = ctc_threshold(1e18, 1e3, 1.0)["margin"]
    m2 = ctc_threshold(1e18, 1e3, 10.0)["margin"]
    assert m2 > m1 * 100  # R^3 scaling
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt007_b2.py -v")
    print("  python -m zt005.run_zt007_b2 --mode quick")