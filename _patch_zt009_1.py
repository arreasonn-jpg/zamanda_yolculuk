"""Auto-patch for ZT-009.1 — Sensitivity analysis."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt009_1.py"
TST = ROOT / "tests" / "test_zt009_1.py"

SRC_CODE = r'''
"""ZT-009.1 - Sensitivity analysis on h_00 across device parameters."""
import argparse, json
from pathlib import Path
import numpy as np

# Scaling laws:
#   h_00 ~ I^2 * f^a * R^b * eta^c * geometry_factor / D
# From our numerical runs:
#   - I^2 (confirmed ZT-006.4)
#   - R^2 (roughly, from B ~ mu0 I R^2 / ...)
#   - f (weak, through omega factor, ~f^0 with quasi-static)
#   - eta (linear in efficiency)
#   - D (inverse from Poisson solver dimensions)

I_BASE = 100.0        # A
R_BASE = 0.12         # m (coil radius)
F_BASE = 1.0e6        # Hz
ETA_BASE = 0.85       # efficiency
D_BASE = 2.0          # m (temporal region diameter)
H00_BASE = 5.9e-48    # reference value at G1, 100 A


def h00_model(I=I_BASE, R=R_BASE, f=F_BASE, eta=ETA_BASE, D=D_BASE):
    """
    Empirical scaling model calibrated to ZT-006.3.1 baseline.
    """
    return H00_BASE * (I / I_BASE)**2 * (R / R_BASE)**2 * \
           (f / F_BASE)**0.0 * (eta / ETA_BASE)**1.0 * (D_BASE / D)


def sensitivity_index(param_name, rel_delta=0.10, **kwargs):
    """
    Compute normalized sensitivity: (dh/h) / (dp/p).
    """
    h0 = h00_model(**kwargs)
    p0 = kwargs.get(param_name, None)
    if p0 is None:
        # get default
        defaults = {"I": I_BASE, "R": R_BASE, "f": F_BASE,
                    "eta": ETA_BASE, "D": D_BASE}
        p0 = defaults[param_name]

    new_kwargs = dict(kwargs)
    new_kwargs[param_name] = p0 * (1 + rel_delta)
    h1 = h00_model(**new_kwargs)
    return ((h1 - h0) / h0) / rel_delta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 84)
    print("ZAMANDA YOLCULUK - ZT-009.1")
    print("SENSITIVITY ANALYSIS")
    print("=" * 84)
    print()

    # --- Sensitivity indices ---
    print("=" * 84)
    print("SENSITIVITY INDICES (normalized)")
    print("=" * 84)
    print("S_i = (dh/h) / (dp/p)  [10% change in p -> ? change in h]")
    print()
    print(f"{'Parameter':<25s} {'Base':<15s} {'S_i':>10s} {'Effect of +10%':>16s}")
    print("-" * 84)

    params = [
        ("I", "Current (A)", I_BASE),
        ("R", "Coil radius (m)", R_BASE),
        ("f", "Frequency (Hz)", F_BASE),
        ("eta", "Efficiency", ETA_BASE),
        ("D", "Region diameter (m)", D_BASE),
    ]
    results = {}
    for key, name, base in params:
        S = sensitivity_index(key, rel_delta=0.10)
        change_pct = S * 10.0
        results[key] = {"S": float(S), "name": name, "base": float(base)}
        print(f"{name:<25s} {base:<15.4e} {S:>10.3f} "
              f"{change_pct:>14.2f}%")

    print()

    # --- Rank parameters ---
    print("=" * 84)
    print("PARAMETER RANKING (by |S_i|)")
    print("=" * 84)
    ranked = sorted(results.items(), key=lambda x: -abs(x[1]["S"]))
    for i, (key, r) in enumerate(ranked, 1):
        print(f"  {i}. {r['name']:<25s} |S| = {abs(r['S']):.3f}")
    print()

    # --- Best single-parameter intervention ---
    print("=" * 84)
    print("SINGLE-PARAMETER IMPROVEMENT (up to factor of 10)")
    print("=" * 84)
    h0 = h00_model()
    print(f"Baseline h_00 = {h0:.4e}")
    print()
    print(f"{'Parameter':<25s} {'10x':>15s} {'100x':>15s} {'10^6x':>15s}")
    print("-" * 84)
    for key, r in ranked:
        row = []
        for scale in [10, 100, 1e6]:
            new_kwargs = {}
            defaults = {"I": I_BASE, "R": R_BASE, "f": F_BASE,
                        "eta": ETA_BASE, "D": D_BASE}
            for k, v in defaults.items():
                new_kwargs[k] = v * (scale if k == key else 1.0)
            # eta can't exceed 1.0
            if key == "eta":
                new_kwargs["eta"] = 1.0
            h_new = h00_model(**new_kwargs)
            row.append(f"{h_new:.2e}")
        print(f"{r['name']:<25s} {row[0]:>15s} {row[1]:>15s} {row[2]:>15s}")
    print()

    # --- Gap to target ---
    print("=" * 84)
    print("HEDEF: 1 s / yil  =>  h_00_target ~ 6.3e-8")
    print("=" * 84)
    h_target = 2.0 / 3.15576e7
    gap = h_target / h0
    print(f"  h_00 baseline:  {h0:.4e}")
    print(f"  h_00 target:    {h_target:.4e}")
    print(f"  Gap:            10^{np.log10(gap):.2f}")
    print()
    print("  Tek parametre ile (10^6x scale):")
    for key, r in ranked:
        defaults = {"I": I_BASE, "R": R_BASE, "f": F_BASE,
                    "eta": ETA_BASE, "D": D_BASE}
        new_kwargs = dict(defaults)
        new_kwargs[key] = defaults[key] * 1e6
        if key == "eta":
            new_kwargs[key] = 1.0
        h_new = h00_model(**new_kwargs)
        g = h_target / h_new
        print(f"    {r['name']:<25s} h_new = {h_new:.3e}  "
              f"gap = 10^{np.log10(g):.2f}")
    print()

    # --- Conclusion ---
    print("=" * 84)
    print("SONUC")
    print("=" * 84)
    top = ranked[0]
    print(f"  En etkili parametre: {top[1]['name']}")
    print(f"  |S| = {abs(top[1]['S']):.3f}")
    print()
    print("  Bu, h_00'in hangi fiziksel buyukluge en duyarli oldugunu")
    print("  gosterir. Ancak 10^40 gap'i kapatmak icin TEK parametre")
    print("  yetmez - tum parametrelerin es zamanli optimize edilmesi")
    print("  ve muhtemelen yeni fizik gerekir.")
    print()

    out = Path("zt009_1_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-009.1",
        "mode": a.mode,
        "h00_baseline": float(h0),
        "h00_target_1s_per_year": float(h_target),
        "gap_log10": float(np.log10(gap)),
        "sensitivity": results,
        "ranking": [k for k, _ in ranked],
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: sensitivity analysis; no new physics.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-009.1 tests."""
import numpy as np
from zt005.run_zt009_1 import (
    h00_model, sensitivity_index,
    I_BASE, R_BASE, H00_BASE,
)


def test_baseline_recovers_H00():
    h = h00_model()
    assert abs(h - H00_BASE) / H00_BASE < 1e-12


def test_current_quadratic():
    h1 = h00_model(I=I_BASE)
    h2 = h00_model(I=2*I_BASE)
    assert abs(h2 / h1 - 4.0) < 1e-12


def test_radius_quadratic():
    h1 = h00_model(R=R_BASE)
    h2 = h00_model(R=2*R_BASE)
    assert abs(h2 / h1 - 4.0) < 1e-12


def test_diameter_inverse():
    h1 = h00_model(D=1.0)
    h2 = h00_model(D=2.0)
    assert abs(h1 / h2 - 2.0) < 1e-12


def test_sensitivity_current_is_2():
    """S_I = dh/dI * I/h = 2 for quadratic."""
    S = sensitivity_index("I", rel_delta=0.001)
    assert abs(S - 2.0) < 1e-3


def test_sensitivity_radius_is_2():
    S = sensitivity_index("R", rel_delta=0.001)
    assert abs(S - 2.0) < 1e-3


def test_sensitivity_diameter_is_minus_1():
    S = sensitivity_index("D", rel_delta=0.001)
    assert abs(S - (-1.0)) < 1e-3


def test_sensitivity_frequency_is_zero():
    """h does not depend on f in our model."""
    S = sensitivity_index("f", rel_delta=0.001)
    assert abs(S) < 1e-3
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt009_1.py -v")
    print("  python -m zt005.run_zt009_1 --mode quick")