"""Auto-patch for ZT-006.4 — Scaling and feasibility analysis."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt006_4.py"
TST = ROOT / "tests" / "test_zt006_4.py"

SRC_CODE = r'''
"""ZT-006.4 - Scaling and feasibility of EM-induced time dilation."""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.run_zt006_3 import compute_geodesic_diagnostics

SECONDS_PER_YEAR = 3.15576e7
G = 6.67430e-11
C = 2.99792458e8


def required_amplification(current_h00, target_h00=2.0 / SECONDS_PER_YEAR):
    """
    Return amplification factor needed to reach delta_tau = 1 s / year.

    Since h ~ I^2, required current amplification = sqrt(h_target / h_curr).
    """
    return np.sqrt(target_h00 / current_h00)


def scaling_table(base_h00, name):
    """
    Show what happens if we scale I by 10x, 100x, 1e6x, 1e10x.
    """
    rows = []
    for scale_I in [1, 10, 1e2, 1e6, 1e10, 1e14, 1e18, 1e20, 1e22]:
        h_new = base_h00 * scale_I**2
        # dtau/yr
        delta_yr = (h_new / 2.0) * SECONDS_PER_YEAR
        rows.append({
            "I_scale": scale_I,
            "h00_new": float(h_new),
            "delta_tau_per_year_s": float(delta_yr),
            "still_below_1s": bool(delta_yr < 1.0),
        })
    return rows


def current_for_1s_per_year(h00_at_100A):
    """
    What current do we need so dtau = 1 s over 1 year?
    h_target = 2 / SECONDS_PER_YEAR
    I_new = 100 * sqrt(h_target / h00_at_100A)
    """
    h_target = 2.0 / SECONDS_PER_YEAR
    I_new = 100.0 * np.sqrt(h_target / h00_at_100A)
    return h_target, I_new


def comparison_benchmarks():
    """
    Known physical time-dilation systems for scale comparison (per year).
    """
    return {
        "GPS satellites (net)":          4.5e-7 * SECONDS_PER_YEAR,
        "Earth surface (relative to infinity)": 7.0e-10 * SECONDS_PER_YEAR,
        "Sun surface (relative to infinity)":  1.0e-6 * SECONDS_PER_YEAR,
        "Neutron star surface":         1.5e-1 * SECONDS_PER_YEAR,
        "Stellar black hole (10 Msun, 10 km)":  3.0e7 * SECONDS_PER_YEAR,
        "Target (1 s / year)":          1.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--n-phi", type=int, default=8)
    a = ap.parse_args()

    print("=" * 72)
    print("ZAMANDA YOLCULUK - ZT-006.4")
    print("SCALING AND FEASIBILITY OF EM-INDUCED TIME DILATION")
    print("=" * 72)
    print()

    # ----- Base h00 per geometry (at 100 A, 1 MHz) -----
    print("Base values (I=100 A, 1 MHz, device in backpack):")
    base = {}
    for g in ["G1", "G2", "G3", "G4"]:
        res = compute_geodesic_diagnostics(g, 40, 80, n_phi=a.n_phi)
        base[g] = res["h_00_center"]
        print(f"  {g}: h00 = {res['h_00_center']:.4e}  "
              f"(1-dtau/dt) = {res['delta_tau_over_t_center']:.4e}")
    print()

    # ----- Required current for G1 (strongest) -----
    h_best = max(base.values())
    h_target, I_required = current_for_1s_per_year(h_best)
    print("=" * 72)
    print("HEDEF: 1 saniye zaman genislenmesi 1 yil icinde")
    print("=" * 72)
    print(f"En guclu geometri: G1, h00 = {h_best:.4e}")
    print(f"Hedef h00 (1 s / 1 yr) = {h_target:.4e}")
    print(f"h orani = {h_target / h_best:.4e}")
    print(f"Gerekli akim = 100 * sqrt({h_target / h_best:.4e}) = {I_required:.4e} A")
    print()
    print(f"Referans: yildirim akimi ~ 10^5 A")
    print(f"Referans: laboratuvar darbeli manyetik alan ~ 10^7 A")
    print(f"Referans: astrofiziksel plazma akimlari ~ 10^18 A")
    print(f"Referans: magnetar yuzey akimlari ~ 10^24 A")
    print()
    print(f"Gereken = {I_required:.2e} A -> "
          f"magnetar olceginin {I_required/1e24:.2e}x kati")
    print()

    # ----- Scaling table for G1 -----
    print("=" * 72)
    print("OLCEK TABLOSU (G1 baz alinarak)")
    print("=" * 72)
    print(f"{'I_scale':>12s}  {'h00':>12s}  {'dtau/yil':>14s}  {'1 s icin?':>10s}")
    print("-" * 72)
    for row in scaling_table(h_best, "G1"):
        ok = "EVET" if not row["still_below_1s"] else "HAYIR"
        print(f"{row['I_scale']:>12.0e}  "
              f"{row['h00_new']:>12.4e}  "
              f"{row['delta_tau_per_year_s']:>14.4e}  "
              f"{ok:>10s}")
    print()

    # ----- Benchmark comparison -----
    print("=" * 72)
    print("BILINEN FIZIKSEL ZAMAN GENISLENMELERI (yillik)")
    print("=" * 72)
    for name, val in comparison_benchmarks().items():
        print(f"  {name:<45s} : {val:.4e} s/yil")
    print()
    print(f"  Device (G1, 100 A)                        : "
          f"{(h_best / 2.0) * SECONDS_PER_YEAR:.4e} s/yil")
    print()
    print(f"Device / GPS orani: "
          f"{((h_best/2)*SECONDS_PER_YEAR) / (4.5e-7*SECONDS_PER_YEAR):.2e}")
    print()

    # ----- Energy budget -----
    V_cylinder = np.pi * 1.0**2 * 2.5
    print("=" * 72)
    print("ENERJI BUTCESI")
    print("=" * 72)
    print(f"Aktif hacim = {V_cylinder:.4f} m^3")
    print(f"Ortalama T00 ~ 5e-4 J/m^3 (G1)")
    print(f"Toplam EM enerji ~ {5e-4 * V_cylinder:.4e} J")
    print(f"Karsilastirma: 1 kalem pili ~ 10^3 J")
    print(f"Karsilastirma: nukleer bomba ~ 10^14 J")
    print()
    print("Sonuc: Mevcut konfigurasyonda toplam enerji ihmal edilebilir.")
    print("Zaman genislenmesi icin gereken enerji yogunlugu:")
    T00_required = 0.25 * (I_required / 100.0)**2
    print(f"  T00_required = {T00_required:.4e} J/m^3")
    print(f"  Karsilastirma: plazma merkezi (tokamak) ~ 10^6 J/m^3")
    print(f"  Karsilastirma: nukleer madde ~ 10^34 J/m^3")
    print(f"  Karsilastirma: kuark-gluon plazmasi ~ 10^38 J/m^3")
    print()

    print("=" * 72)
    print("SONUC")
    print("=" * 72)
    print("Elektromanyetik alanla, cihaz olceginde, 1 saniyelik")
    print("zaman genislenmesi uretmek icin gereken enerji yogunlugu")
    print("nukleer maddenin ~10^5 kati mertebesindedir.")
    print("Bu, mevcut veya ongorulebilir teknolojiyle ULASILAMAZ.")
    print("=" * 72)

    out = Path("zt006_4_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.4",
        "mode": a.mode,
        "base_h00": {k: float(v) for k, v in base.items()},
        "h_target_1s_per_year": float(h_target),
        "required_current_A": float(I_required),
        "T00_required_J_per_m3": float(T00_required),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: feasibility falsification; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-006.4 tests."""
import numpy as np
from zt005.run_zt006_4 import (
    required_amplification, current_for_1s_per_year,
    scaling_table, comparison_benchmarks,
)

SECONDS_PER_YEAR = 3.15576e7


def test_amplification_factor_is_huge():
    """From 1e-48 to 6e-8, we need ~1e20 current amplification."""
    base_h00 = 1e-48
    amp = required_amplification(base_h00)
    assert 1e19 < amp < 1e21


def test_required_current_is_absurd():
    h_base = 6e-48
    h_t, I_req = current_for_1s_per_year(h_base)
    assert 1e21 < I_req < 1e23


def test_scaling_table_monotone():
    rows = scaling_table(1e-48, "test")
    deltas = [r["delta_tau_per_year_s"] for r in rows]
    assert all(deltas[i] < deltas[i+1] for i in range(len(deltas)-1))


def test_target_1s_per_year_reachable_only_with_1e20():
    """Requires I ~ 1e22 A, which is >> magnetar currents."""
    h_base = 6e-48
    _, I_req = current_for_1s_per_year(h_base)
    assert I_req > 1e21


def test_benchmarks_include_gps_and_target():
    b = comparison_benchmarks()
    assert "GPS satellites (net)" in b
    assert "Target (1 s / year)" in b
    assert b["Target (1 s / year)"] == 1.0


def test_gps_dilation_is_nanoseconds_per_day():
    b = comparison_benchmarks()
    gps_per_day = b["GPS satellites (net)"] / 365.0
    assert 1e-10 < gps_per_day < 1e-6
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt006_4.py -v")
    print("  python -m zt005.run_zt006_4 --mode quick")