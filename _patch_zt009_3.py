"""Auto-patch for ZT-009.3 — Best-case engineering constraints."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt009_3.py"
TST = ROOT / "tests" / "test_zt009_3.py"

SRC_CODE = r'''
"""ZT-009.3 - Best-case h_00 under engineering constraints.

Push every parameter to its physical/engineering limit and see
how close to the target we can get.
"""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.run_zt009_1 import h00_model, H00_BASE, \
    I_BASE, R_BASE, ETA_BASE, D_BASE

SECONDS_PER_YEAR = 3.15576e7
H00_TARGET = 2.0 / SECONDS_PER_YEAR


# --- Engineering constraints ---
ENVELOPE = {
    "max_height_m": 0.45,
    "max_width_m": 0.32,
    "max_depth_m": 0.15,
    "max_coil_radius_m": 0.10,  # fits with margin
}

THERMAL = {
    "max_temp_C": 40,
    "ambient_C": 25,
    "copper_density_A_per_mm2": 10,      # conservative
    "cooled_copper_A_per_mm2": 100,      # liquid-cooled
    "superconductor_A_per_mm2": 1000,    # HTS at 77 K
}

POWER = {
    "battery_wh": 200,
    "min_runtime_min": 30,
}

# Stacking: N coils in series along z, each contributes
# h ~ N^2 * (single-coil h) if phase-locked
MAX_STAGES = 10

# Resonance: high-Q cavity (superconducting niobium)
MAX_Q = 1e4


def tier_scenarios():
    """
    Three tiers of engineering feasibility.
    """
    return [
        {
            "name": "Baseline (copper, 100 A)",
            "I_A": 100,
            "R_m": 0.07,
            "eta": 0.85,
            "N_stages": 1,
            "Q_factor": 1,
            "cooling": "passive",
            "cost": "low",
            "tr_level": 6,
        },
        {
            "name": "Advanced (cooled copper, 500 A)",
            "I_A": 500,
            "R_m": 0.09,
            "eta": 0.92,
            "N_stages": 4,
            "Q_factor": 100,
            "cooling": "liquid",
            "cost": "medium",
            "tr_level": 6,
        },
        {
            "name": "Best-case (superconducting HTS)",
            "I_A": 10000,
            "R_m": 0.10,
            "eta": 0.99,
            "N_stages": MAX_STAGES,
            "Q_factor": MAX_Q,
            "cooling": "cryogenic",
            "cost": "high",
            "tr_level": 8,
        },
    ]


def compute_h00_with_engineering(
    I_A, R_m, eta, N_stages=1, Q_factor=1
):
    """
    h_00 = base_model(I, R, eta) * N_stages^2 * Q_factor
    """
    h_single = h00_model(I=I_A, R=R_m, eta=eta, D=D_BASE)
    return h_single * (N_stages**2) * Q_factor


def max_current_thermal(
    wire_area_mm2=1.0, current_density=None, cooling="passive"
):
    """
    Maximum current from wire current density.
    """
    if current_density is None:
        current_density = THERMAL["copper_density_A_per_mm2"]
        if cooling == "liquid":
            current_density = THERMAL["cooled_copper_A_per_mm2"]
        elif cooling == "cryogenic":
            current_density = THERMAL["superconductor_A_per_mm2"]
    return wire_area_mm2 * current_density


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 96)
    print("ZAMANDA YOLCULUK - ZT-009.3")
    print("BEST-CASE ENGINEERING CONSTRAINTS")
    print("=" * 96)
    print()

    # --- Engineering envelope ---
    print("=" * 96)
    print("ENGINEERING ENVELOPE")
    print("=" * 96)
    for k, v in ENVELOPE.items():
        print(f"  {k:<28s} {v}")
    print()
    for k, v in THERMAL.items():
        print(f"  {k:<28s} {v}")
    print()
    for k, v in POWER.items():
        print(f"  {k:<28s} {v}")
    print()

    # --- Scenarios ---
    print("=" * 96)
    print("TIER SCENARIOS")
    print("=" * 96)
    print(f"{'Name':<38s} {'I (A)':>10s} {'R (m)':>8s} "
          f"{'eta':>6s} {'N':>4s} {'Q':>8s} {'TR':>4s}")
    print("-" * 96)
    for t in tier_scenarios():
        print(f"{t['name']:<38s} {t['I_A']:>10.0f} "
              f"{t['R_m']:>8.3f} {t['eta']:>6.3f} "
              f"{t['N_stages']:>4d} {t['Q_factor']:>8.0f} "
              f"{t['tr_level']:>4d}")
    print()

    # --- Compute h_00 per tier ---
    print("=" * 96)
    print("COMPUTED h_00 PER TIER")
    print("=" * 96)
    print(f"{'Tier':<38s} {'h_00':>14s} {'Rel to base':>14s} "
          f"{'Gap to 1s/yr':>14s}")
    print("-" * 96)
    results = []
    for t in tier_scenarios():
        h = compute_h00_with_engineering(
            t["I_A"], t["R_m"], t["eta"],
            t["N_stages"], t["Q_factor"],
        )
        rel = h / H00_BASE
        gap_log = np.log10(H00_TARGET / h)
        results.append({
            **t,
            "h00": float(h),
            "rel_to_base": float(rel),
            "gap_log10": float(gap_log),
        })
        print(f"{t['name']:<38s} {h:>14.4e} {rel:>14.4e} "
              f"10^{gap_log:>4.1f}")
    print()

    # --- Best case ---
    best = max(results, key=lambda r: r["h00"])
    print("=" * 96)
    print("BEST-CASE SCENARIO")
    print("=" * 96)
    print(f"  Configuration:  {best['name']}")
    print(f"  I (A):          {best['I_A']:.2e}")
    print(f"  R (m):          {best['R_m']:.3f}")
    print(f"  eta:            {best['eta']:.3f}")
    print(f"  Stages:         {best['N_stages']}")
    print(f"  Q factor:       {best['Q_factor']:.0e}")
    print(f"  h_00:           {best['h00']:.4e}")
    print(f"  Gap to target:  10^{best['gap_log10']:.1f}")
    print()

    # --- Max possible current from wire (thermal) ---
    print("=" * 96)
    print("MAXIMUM CURRENT FROM WIRE HEATING (physical bound)")
    print("=" * 96)
    for cooling in ["passive", "liquid", "cryogenic"]:
        for area_mm2 in [0.5, 1.0, 5.0, 20.0]:
            I_max = max_current_thermal(area_mm2, cooling=cooling)
            print(f"  {cooling:<10s} A={area_mm2:>5.1f} mm^2  "
                  f"->  I_max = {I_max:>10.0f} A")
    print()

    # --- Extrapolate: how far would we have to go? ---
    print("=" * 96)
    print("EXTRAPOLATION: how far would we need to go?")
    print("=" * 96)
    h_best = best["h00"]
    gap_linear = H00_TARGET / h_best
    print(f"  Best-case h_00:      {h_best:.4e}")
    print(f"  Target h_00:         {H00_TARGET:.4e}")
    print(f"  Linear gap:          {gap_linear:.4e}")
    print()
    print(f"  If we could scale current only (h ~ I^2):")
    I_needed = best["I_A"] * np.sqrt(gap_linear)
    print(f"    I_needed = {I_needed:.4e} A")
    print(f"    Ref: magnetar current ~ 1e24 A")
    print()
    print(f"  If we could scale radius only (h ~ R^2):")
    R_needed = best["R_m"] * np.sqrt(gap_linear)
    print(f"    R_needed = {R_needed:.4e} m")
    print(f"    Ref: solar system size ~ 1e13 m")
    print()
    print(f"  If we could scale stages only (h ~ N^2):")
    N_needed = best["N_stages"] * np.sqrt(gap_linear)
    print(f"    N_needed = {N_needed:.4e}")
    print()
    print(f"  If we could scale Q only (h ~ Q):")
    Q_needed = best["Q_factor"] * gap_linear
    print(f"    Q_needed = {Q_needed:.4e}")
    print(f"    Ref: best lab cavity Q ~ 1e12")
    print()

    # --- Conclusion ---
    print("=" * 96)
    print("SONUC")
    print("=" * 96)
    print(f"  Even the best-case configuration (superconducting, 10^4 A,")
    print(f"  R = 10 cm, 10 stages, Q = 10^4) gives h_00 = {best['h00']:.2e}")
    print(f"  Still 10^{best['gap_log10']:.1f} short of the 1 s/year target.")
    print()
    print("  No combination of currently conceivable engineering")
    print("  improvements can bridge this gap. New physics is required.")
    print()

    # --- Figure ---
    fig_path = Path("zt009_3_best_case.png")
    fig_ok = save_best_case_figure(results, str(fig_path))
    print(f"Figure saved: {fig_path} ({'OK' if fig_ok else 'FAIL'})")
    print()

    out = Path("zt009_3_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-009.3",
        "mode": a.mode,
        "envelope": ENVELOPE,
        "thermal": THERMAL,
        "power": POWER,
        "tiers": results,
        "best_case": best,
        "extrapolation": {
            "I_needed_A": float(I_needed),
            "R_needed_m": float(R_needed),
            "N_needed": float(N_needed),
            "Q_needed": float(Q_needed),
        },
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: engineering-bound analysis; no new physics.")


def save_best_case_figure(results, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        names = [r["name"].split("(")[0].strip() for r in results]
        log_h = [np.log10(r["h00"]) for r in results]
        log_target = np.log10(H00_TARGET)
        log_base = np.log10(H00_BASE)

        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.barh(names, log_h, color=["#999999",
                                            "#f39c12", "#2ca02c"])
        ax.axvline(log_target, color="red", linestyle="--",
                   label=f"Target (1 s/yr): 10^{log_target:.1f}")
        ax.axvline(log_base, color="blue", linestyle=":",
                   label=f"Baseline: 10^{log_base:.1f}")

        for bar, lg in zip(bars, log_h):
            ax.text(lg - 1, bar.get_y() + bar.get_height()/2,
                    f"10^{lg:.1f}", va="center", ha="right",
                    fontsize=10, color="white")

        ax.set_xlabel("log10(h_00)")
        ax.set_title("Best-case h_00 across engineering tiers vs target")
        ax.set_xlim(log_base - 10, log_target + 5)
        ax.legend(loc="lower right")
        ax.grid(axis="x", alpha=0.3)
        plt.tight_layout()
        fig.savefig(out_path, dpi=120)
        plt.close(fig)
        return True
    except Exception as e:
        print(f"Figure save failed: {e}")
        return False


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-009.3 tests."""
import numpy as np
from zt005.run_zt009_3 import (
    tier_scenarios, compute_h00_with_engineering,
    max_current_thermal, ENVELOPE, THERMAL,
    H00_TARGET, MAX_STAGES, MAX_Q,
)


def test_three_tiers():
    t = tier_scenarios()
    assert len(t) == 3


def test_tiers_increasing_power():
    t = tier_scenarios()
    currents = [x["I_A"] for x in t]
    assert currents[0] < currents[1] < currents[2]


def test_best_case_uses_superconducting():
    t = tier_scenarios()
    best = t[-1]
    assert best["I_A"] >= 1000
    assert best["eta"] > 0.95
    assert best["N_stages"] >= 5


def test_compute_h00_scales_quadratically_in_current():
    h1 = compute_h00_with_engineering(100, 0.1, 0.9)
    h2 = compute_h00_with_engineering(200, 0.1, 0.9)
    assert abs(h2 / h1 - 4.0) < 1e-9


def test_compute_h00_scales_quadratically_in_stages():
    h1 = compute_h00_with_engineering(100, 0.1, 0.9, N_stages=1)
    h2 = compute_h00_with_engineering(100, 0.1, 0.9, N_stages=10)
    assert abs(h2 / h1 - 100.0) < 1e-9


def test_compute_h00_scales_linearly_in_Q():
    h1 = compute_h00_with_engineering(100, 0.1, 0.9, Q_factor=1)
    h2 = compute_h00_with_engineering(100, 0.1, 0.9, Q_factor=1000)
    assert abs(h2 / h1 - 1000.0) < 1e-9


def test_max_current_scales_with_density():
    I1 = max_current_thermal(1.0, cooling="passive")
    I2 = max_current_thermal(1.0, cooling="cryogenic")
    assert I2 > I1 * 10


def test_best_case_still_far_from_target():
    t = tier_scenarios()
    best = t[-1]
    h = compute_h00_with_engineering(
        best["I_A"], best["R_m"], best["eta"],
        best["N_stages"], best["Q_factor"],
    )
    # Even best case must be far below target
    assert h < H00_TARGET / 1e20
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt009_3.py -v")
    print("  python -m zt005.run_zt009_3 --mode quick")