"""Auto-patch for ZT-007B.5 — Alcubierre warp drive feasibility."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt007_b5.py"
TST = ROOT / "tests" / "test_zt007_b5.py"

SRC_CODE = r'''
"""ZT-007B.5 - Alcubierre warp drive: exotic energy requirement."""
import argparse, json
from pathlib import Path
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
HBAR = 1.054571817e-34
MSUN = 1.98892e30
MEARTH = 5.972e24
MJUP = 1.898e27


def alcubierre_energy(v_s, R):
    """
    Total negative energy of an Alcubierre warp bubble (order of magnitude):

        E ~ - v_s^2 R^2 c^4 / (12 G)

    Reference: Alcubierre 1994; Pfenning-Ford 1997.

    Parameters
    ----------
    v_s : float  bubble velocity (m/s)
    R   : float  bubble radius (m)

    Returns
    -------
    E_negative : float  in Joules (negative)
    """
    return -v_s**2 * R**2 * C**4 / (12.0 * G)


def alcubierre_energy_density(v_s, R):
    """
    Typical negative energy density in the bubble wall (peak):

        rho ~ - v_s^2 c^2 / (8 pi G R^2)
    """
    return -v_s**2 * C**2 / (8.0 * np.pi * G * R**2)


def casimir_max_density(d=1e-9):
    """Casimir |rho| at d nm separation."""
    return np.pi**2 * HBAR * C / (720.0 * d**4)


def mass_equivalent(E):
    """Return |E|/c^2 in kg."""
    return abs(E) / C**2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 78)
    print("ZAMANDA YOLCULUK - ZT-007B.5")
    print("ALCUBIERRE WARP DRIVE: EXOTIC ENERGY REQUIREMENT")
    print("=" * 78)
    print()
    print("Formula: E_neg ~ -v_s^2 R^2 c^4 / (12 G)")
    print()

    # --- Energy vs velocity and radius ---
    print(f"{'v_s/c':>8s} {'R (m)':>10s} {'E_neg (J)':>16s} "
          f"{'|E|/c^2 (kg)':>16s} {'Solar masses':>16s}")
    print("-" * 78)
    for beta in [0.01, 0.1, 1.0, 10.0]:
        v_s = beta * C
        for R in [1.0, 10.0, 100.0, 1e3]:
            E = alcubierre_energy(v_s, R)
            m = mass_equivalent(E)
            m_sun = m / MSUN
            print(f"{beta:>8.3f} {R:>10.1f} {E:>16.4e} "
                  f"{m:>16.4e} {m_sun:>16.4e}")
    print()

    # --- Energy density vs radius at v_s = 0.1 c ---
    print("=" * 78)
    print("ENERGY DENSITY AT v_s = 0.1 c")
    print("=" * 78)
    beta = 0.1
    v_s = beta * C
    print(f"{'R (m)':>10s} {'|rho| (J/m^3)':>22s} "
          f"{'Gap vs Casimir(1nm)':>24s}")
    print("-" * 78)
    rho_cas = casimir_max_density(1e-9)
    for R in [1e-3, 1e-2, 0.1, 1.0, 10.0]:
        rho = abs(alcubierre_energy_density(v_s, R))
        gap = rho / rho_cas
        print(f"{R:>10.1e} {rho:>22.4e} {gap:>24.4e}")
    print()
    print(f"Casimir |rho| at 1 nm = {rho_cas:.4e} J/m^3")
    print()

    # --- Comparison to known objects ---
    print("=" * 78)
    print("KARSILASTIRMA: 100 m yaricap, subluminal warp")
    print("=" * 78)
    R = 100.0
    for beta in [0.001, 0.01, 0.1, 1.0]:
        E = alcubierre_energy(beta * C, R)
        m = mass_equivalent(E)
        print(f"  v={beta:.3f}c, R={R} m  ->  |E| = {abs(E):.4e} J  "
              f"= {m/MSUN:.4e} M_sun  = {m/MJUP:.4e} M_Jup")
    print()

    # --- Pfenning-Ford quantum inequality ---
    print("=" * 78)
    print("PFENNING-FORD BOUND (quantum inequality on warp bubble)")
    print("=" * 78)
    print("For Alcubierre warp bubble, Pfenning-Ford (1997) showed:")
    print("  Delta_E ~ -v_s^2 R^2 c^4 / (12 G)  (Alcubierre's estimate)")
    print("  AND quantum inequality requires tau < R / c  (light-crossing)")
    print("  AND total negative energy >= -3 hbar c / (32 pi^2 tau^2 R)")
    print()
    print("Practical consequence: warp bubble must be thinner than ~")
    print("Planck length to satisfy quantum inequality, which is unphysical.")
    print()

    # --- Conclusion ---
    print("=" * 78)
    print("SONUC")
    print("=" * 78)
    print("Alcubierre warp drive exotic energy requirements (order of magnitude):")
    print()
    print(f"  Lab-scale (R=1 m, v=0.1c)  : |E| ~ {abs(alcubierre_energy(0.1*C, 1.0)):.2e} J")
    print(f"                             :       ~ {mass_equivalent(alcubierre_energy(0.1*C, 1.0))/MEARTH:.2e} Earth masses")
    print()
    print(f"  Subluminal (R=100 m, v=0.1c): |E| ~ {abs(alcubierre_energy(0.1*C, 100.0)):.2e} J")
    print(f"                              :       ~ {mass_equivalent(alcubierre_energy(0.1*C, 100.0))/MJUP:.2e} Jupiter masses")
    print()
    print(f"  Superluminal (R=100 m, v=c) : |E| ~ {abs(alcubierre_energy(1.0*C, 100.0)):.2e} J")
    print(f"                              :       ~ {mass_equivalent(alcubierre_energy(1.0*C, 100.0))/MSUN:.2e} Solar masses")
    print()
    print("Bu, device-scale'de UYGULANAMAZ. Pfenning-Ford quantum")
    print("inequality, negatif enerjiyi kisa sureli sinirlayarak")
    print("bubble'in stabilitesini daha da zorlastirir.")
    print("=" * 78)

    out = Path("zt007_b5_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-007B.5",
        "mode": a.mode,
        "alcubierre_coefficient": C**4 / (12 * G),
        "casimir_max_1nm": float(casimir_max_density(1e-9)),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: Alcubierre feasibility; no warp construction.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-007B.5 tests - Alcubierre."""
import numpy as np
from zt005.run_zt007_b5 import (
    alcubierre_energy, alcubierre_energy_density,
    casimir_max_density, mass_equivalent,
    C, MSUN,
)


def test_energy_negative():
    E = alcubierre_energy(0.1 * C, 10.0)
    assert E < 0


def test_energy_scales_v_squared():
    E1 = abs(alcubierre_energy(0.1 * C, 10.0))
    E2 = abs(alcubierre_energy(0.2 * C, 10.0))
    assert abs(E2 / E1 - 4.0) < 1e-10


def test_energy_scales_R_squared():
    E1 = abs(alcubierre_energy(0.1 * C, 10.0))
    E2 = abs(alcubierre_energy(0.1 * C, 20.0))
    assert abs(E2 / E1 - 4.0) < 1e-10


def test_1m_bubble_requires_earth_mass():
    # 1 m bubble at 0.1c: |E| ~ 1e44 J = ~ 1000 Earth masses
    E = abs(alcubierre_energy(0.1 * C, 1.0))
    m = mass_equivalent(E)
    m_earth = m / 5.972e24
    assert 1e2 < m_earth < 1e4


def test_100m_superluminal_requires_solar_mass():
    # 100 m bubble at c: |E| ~ 1e48 J = ~ 10 Solar masses
    E = abs(alcubierre_energy(1.0 * C, 100.0))
    m = mass_equivalent(E)
    m_sun = m / MSUN
    assert 1.0 < m_sun < 100.0


def test_energy_density_negative():
    rho = alcubierre_energy_density(0.1 * C, 10.0)
    assert rho < 0


def test_casimir_gap_huge():
    """Warp bubble needs ~1e35+ more negative energy than Casimir."""
    rho_warp = abs(alcubierre_energy_density(0.1 * C, 1.0))
    rho_cas  = casimir_max_density(1e-9)
    gap = rho_warp / rho_cas
    assert gap > 1e20
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt007_b5.py -v")
    print("  python -m zt005.run_zt007_b5 --mode quick")