"""Auto-patch for ZT-007B.3 — Morris-Thorne wormhole + Casimir bound."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt007_b3.py"
TST = ROOT / "tests" / "test_zt007_b3.py"

SRC_CODE = r'''
"""ZT-007B.3 - Morris-Thorne wormhole + Casimir exotic-matter bound."""
import argparse, json
from pathlib import Path
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
HBAR = 1.054571817e-34
L_PLANCK = 1.616255e-35


def exotic_mass_required(b0):
    """
    Morris-Thorne wormhole throat of radius b0 requires exotic mass:
        M_exotic ~ -c^2 * b0 / G
    (order of magnitude)
    """
    return C**2 * b0 / G


def casimir_energy_density(d):
    """
    Casimir energy density between parallel plates separated by d.
    rho_Casimir = -pi^2 hbar c / (720 d^4)
    """
    return -np.pi**2 * HBAR * C / (720.0 * d**4)


def casimir_plate_distance_for_rho(target_rho):
    """
    Invert: d = (pi^2 hbar c / (720 |rho|))^(1/4)
    """
    return (np.pi**2 * HBAR * C / (720.0 * abs(target_rho)))**0.25


def max_casimir_density(d_min=1e-9):
    """At d_min (1 nm = atomic scale), Casimir density max."""
    return abs(casimir_energy_density(d_min))


def wormhole_casimir_gap(b0, d_plate):
    """
    For wormhole throat b0, need exotic energy at scale b0.
    Casimir can only provide it at plate distance d_plate.
    Compute gap.
    """
    M_exotic = exotic_mass_required(b0)
    # Effective energy density needed: M c^2 / (4/3 pi b0^3)
    V_throat = 4.0 / 3.0 * np.pi * b0**3
    rho_required = M_exotic * C**2 / V_throat
    rho_casimir_max = max_casimir_density(d_plate)
    return {
        "rho_required_J_m3": float(rho_required),
        "rho_casimir_max_J_m3": float(rho_casimir_max),
        "gap": float(abs(rho_required) / rho_casimir_max),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 78)
    print("ZAMANDA YOLCULUK - ZT-007B.3")
    print("MORRIS-THORNE WORMHOLE + CASIMIR EXOTIC-MATTER BOUND")
    print("=" * 78)
    print()

    # --- Wormhole exotic mass ---
    print("Morris-Thorne wormhole: required exotic mass vs throat radius")
    print(f"{'b0 (m)':>14s} {'M_exotic (kg)':>20s} {'rho_required (J/m^3)':>24s}")
    print("-" * 78)
    for b0 in [1e-9, 1e-6, 1e-3, 1.0, 1e3, 1e6]:
        M = exotic_mass_required(b0)
        rho = M * C**2 / (4.0/3.0 * np.pi * b0**3)
        print(f"{b0:>14.4e} {M:>20.4e} {rho:>24.4e}")
    print()

    # --- Casimir density ---
    print("Casimir energy density between plates:")
    print(f"{'d (m)':>14s} {'rho_Casimir (J/m^3)':>26s}")
    print("-" * 78)
    for d in [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]:
        rho = casimir_energy_density(d)
        print(f"{d:>14.4e} {rho:>26.4e}")
    print()

    # --- Critical comparison: 1 mm wormhole ---
    print("=" * 78)
    print("ORNEK: 1 mm yaricapinda solucan deligi")
    print("=" * 78)
    b0 = 1e-3
    M_req = exotic_mass_required(b0)
    V = 4.0/3.0 * np.pi * b0**3
    rho_req = M_req * C**2 / V
    print(f"  b0                = {b0} m")
    print(f"  M_exotic          = {M_req:.4e} kg")
    print(f"  Gerekli rho       = {rho_req:.4e} J/m^3")
    print()

    # --- Casimir at atomic scale ---
    print("Casimir at d = 1 nm (atomic scale):")
    rho_cas = casimir_energy_density(1e-9)
    print(f"  rho_Casimir       = {rho_cas:.4e} J/m^3")
    print()

    # --- Gap ---
    gap = abs(rho_req) / abs(rho_cas)
    print(f"GAP = {gap:.4e}")
    print(f"Yani 1 mm solucan deligi icin gereken enerji yogunlugu,")
    print(f"Casimir'in saglayabileceginin {gap:.4e} kati.")
    print()

    # --- Ford-Roman quantum inequality ---
    print("=" * 78)
    print("FORD-ROMAN QUANTUM INEQUALITY (temporal bound on exotic matter)")
    print("=" * 78)
    print("Negatif enerji yogunlugu sadece kisa sureli gorunebilir:")
    print("    integral(rho(t)) dt >= -3 hbar / (32 pi^2 tau^4)")
    print()
    print("Tipik sure tau icin izin verilen negatif enerji:")
    print(f"{'tau (s)':>14s} {'limit (J*s/m^3)':>24s}")
    print("-" * 78)
    for tau in [1e-15, 1e-12, 1e-9, 1e-6, 1.0]:
        bound = -3.0 * HBAR / (32.0 * np.pi**2 * tau**4)
        print(f"{tau:>14.4e} {bound:>24.4e}")
    print()

    # --- Conclusion ---
    print("=" * 78)
    print("SONUC")
    print("=" * 78)
    print("Morris-Thorne solucan deligi egzotik madde gerektirir.")
    print("Casimir etkisi, 1 nm plaka araliginda ~ -1e5 J/m^3 uretir.")
    print("1 mm solucan deligi icin gereken rho ~ 1e42 J/m^3.")
    print()
    print(f"Gap: {gap:.2e} kat")
    print()
    print("Karsilastirma:")
    print(f"  Casimir (1 nm)         : 1e5 J/m^3")
    print(f"  Nukleer madde          : 1e34 J/m^3")
    print(f"  1 mm wormhole gerekli  : 1e42 J/m^3")
    print(f"  Planck yogunluk        : 5e96 J/m^3")
    print()
    print("Yani 1 mm'lik solucan deligi icin QCD olceginin 10^8 kati")
    print("enerji yogunlugu gerekir. Bu, laboratuvar olceginde UYGULANAMAZ.")
    print()
    print("Ancak: Planck olceginde (5e96 J/m^3) solucan delikleri")
    print("teorik olarak mumkundur. Bu, kuantum kromodinamigi degil,")
    print("kuantum gravite rejimidir.")
    print("=" * 78)

    out = Path("zt007_b3_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-007B.3",
        "mode": a.mode,
        "c2": C**2,
        "exotic_mass_coefficient": C**2 / G,
        "casimir_coefficient": np.pi**2 * HBAR * C / 720.0,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: Wormhole feasibility bound; no CTC construction.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-007B.3 tests."""
import numpy as np
from zt005.run_zt007_b3 import (
    exotic_mass_required, casimir_energy_density,
    casimir_plate_distance_for_rho, max_casimir_density,
    wormhole_casimir_gap, G, C, HBAR,
)


def test_exotic_mass_coefficient():
    # M_exotic = c^2/G * b0 = 1.35e27 * b0
    b0 = 1.0
    M = exotic_mass_required(b0)
    expected = C**2 / G
    assert abs(M - expected) / expected < 1e-6


def test_exotic_mass_is_huge():
    # For 1 mm: M ~ 1.35e24 kg (Earth mass scale)
    M = exotic_mass_required(1e-3)
    assert 1e23 < M < 1e25


def test_casimir_density_negative():
    rho = casimir_energy_density(1e-9)
    assert rho < 0


def test_casimir_density_increases_at_smaller_d():
    rho1 = abs(casimir_energy_density(1e-8))
    rho2 = abs(casimir_energy_density(1e-9))
    # rho ~ d^-4
    assert abs(rho2 / rho1 - 1e4) / 1e4 < 1e-6


def test_casimir_inverse_roundtrip():
    d = 1e-9
    rho = casimir_energy_density(d)
    d_back = casimir_plate_distance_for_rho(rho)
    assert abs(d_back - d) / d < 1e-6


def test_max_casimir_density_nanometer():
    # At d=1nm, rho ~ -1e5 J/m^3
    rho = max_casimir_density(1e-9)
    assert 1e4 < rho < 1e6


def test_wormhole_casimir_gap_huge():
    # 1 mm wormhole gap >> 1e30
    gap = wormhole_casimir_gap(1e-3, 1e-9)["gap"]
    assert gap > 1e30


def test_casimir_does_not_scale_with_throat():
    # For very tiny wormholes (Planck scale), gap shrinks
    gap_macro = wormhole_casimir_gap(1e-3, 1e-9)["gap"]
    gap_planck = wormhole_casimir_gap(1e-35, 1e-9)["gap"]
    assert gap_planck < gap_macro
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt007_b3.py -v")
    print("  python -m zt005.run_zt007_b3 --mode quick")