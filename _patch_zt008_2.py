"""Auto-patch for ZT-008.2 — Power + thermal budget."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt008_2.py"
TST = ROOT / "tests" / "test_zt008_2.py"

SRC_CODE = r'''
"""ZT-008.2 - Power and thermal budget for backpack."""
import argparse, json
from pathlib import Path
import numpy as np

# Electrical constants
I_OPERATING = 100.0        # A (baseline)
F_DRIVE = 1.0e6            # Hz (1 MHz)
COIL_R_PER_LOOP = 0.005    # ohm (typical copper, small coil)
N_LOOPS = 6
HV_EFFICIENCY = 0.85       # driver + tank efficiency
BATTERY_CAPACITY_WH = 200  # Wh (typical laptop-class)
COOLING_TARGET_T = 40.0    # Celsius
AMBIENT_T = 25.0           # Celsius


def coil_resistance():
    """Total resistance of all loops."""
    return COIL_R_PER_LOOP * N_LOOPS


def joule_heating():
    """
    Joule heating: P = I^2 R (RMS, since AC).
    I here is peak; RMS = I / sqrt(2).
    """
    R = coil_resistance()
    I_rms = I_OPERATING / np.sqrt(2.0)
    return I_rms**2 * R


def power_breakdown():
    """Power consumption per component."""
    p_coil = joule_heating()
    p_driver_overhead = p_coil / HV_EFFICIENCY - p_coil
    p_control = 3.0    # W (MCU + sensors + SBC duty cycle)
    p_safety = 1.0     # W (relay + monitors)
    return {
        "coil_joule": float(p_coil),
        "driver_loss": float(p_driver_overhead),
        "control": p_control,
        "safety": p_safety,
    }


def total_power():
    b = power_breakdown()
    return sum(b.values())


def battery_runtime(capacity_wh=None, power_w=None):
    """Runtime in minutes given capacity and power draw."""
    if capacity_wh is None:
        capacity_wh = BATTERY_CAPACITY_WH
    if power_w is None:
        power_w = total_power()
    return capacity_wh / power_w * 60.0


def thermal_steady_state(P_total, cooling_coeff=1.0):
    """
    Simple thermal model: P = h * A * (T - T_ambient)
    where h*A ~ cooling_coeff (W/K).
    """
    return AMBIENT_T + P_total / cooling_coeff


def required_cooling_coeff(P_total, T_max):
    """Cooling coefficient needed to keep T_max."""
    return P_total / (T_max - AMBIENT_T)


def heat_flux_density(P_total, area_m2=0.45 * 0.32):
    """Heat flux per unit area (W/m^2)."""
    return P_total / area_m2


def save_thermal_figure(P_total, out_path):
    """Plot T vs cooling coefficient."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        hA = np.linspace(0.1, 20.0, 200)
        T_steady = AMBIENT_T + P_total / hA

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(hA, T_steady, "b-", linewidth=2)
        ax.axhline(COOLING_TARGET_T, color="r", linestyle="--",
                   label=f"Target T < {COOLING_TARGET_T} C")
        ax.axhline(80.0, color="darkred", linestyle="--",
                   label="Danger T > 80 C")
        ax.set_xlabel("Cooling coefficient h*A (W/K)")
        ax.set_ylabel("Steady-state temperature (C)")
        ax.set_title(f"Thermal steady-state (P = {P_total:.2f} W)")
        ax.set_ylim(20, 120)
        ax.legend()
        ax.grid(alpha=0.3)
        plt.tight_layout()
        fig.savefig(out_path, dpi=120)
        plt.close(fig)
        return True
    except Exception as e:
        print(f"Figure save failed: {e}")
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 84)
    print("ZAMANDA YOLCULUK - ZT-008.2")
    print("POWER AND THERMAL BUDGET")
    print("=" * 84)
    print()
    print(f"Operating current: {I_OPERATING} A")
    print(f"Drive frequency:   {F_DRIVE/1e6:.1f} MHz")
    print(f"Coil resistance:   {coil_resistance()} ohm (total)")
    print()

    # --- Power breakdown ---
    print("=" * 84)
    print("POWER BREAKDOWN")
    print("=" * 84)
    b = power_breakdown()
    P_tot = sum(b.values())
    for name, val in b.items():
        pct = val / P_tot * 100
        print(f"  {name:<20s} {val:>10.3f} W   ({pct:>5.1f}%)")
    print("-" * 60)
    print(f"  {'TOTAL':<20s} {P_tot:>10.3f} W")
    print()

    # --- Battery runtime ---
    print("=" * 84)
    print("BATTERY RUNTIME")
    print("=" * 84)
    runtime_min = battery_runtime()
    print(f"  Battery capacity: {BATTERY_CAPACITY_WH} Wh")
    print(f"  Total power:      {P_tot:.2f} W")
    print(f"  Runtime:          {runtime_min:.2f} min "
          f"({runtime_min/60:.2f} hours)")
    print()

    # --- Thermal ---
    print("=" * 84)
    print("THERMAL ANALYSIS")
    print("=" * 84)
    T_steady_default = thermal_steady_state(P_tot, cooling_coeff=1.0)
    hA_needed = required_cooling_coeff(P_tot, COOLING_TARGET_T)
    q_density = heat_flux_density(P_tot)

    print(f"  Ambient:           {AMBIENT_T} C")
    print(f"  Target max:        {COOLING_TARGET_T} C")
    print(f"  Total heat load:   {P_tot:.2f} W")
    print(f"  Heat flux:         {q_density:.2f} W/m^2")
    print()
    print(f"  With h*A = 1.0 W/K: T_steady = {T_steady_default:.2f} C")
    print(f"  Required h*A:       {hA_needed:.3f} W/K")
    print()

    # --- Assessment ---
    print("=" * 84)
    print("DEGERLENDIRME")
    print("=" * 84)
    if runtime_min < 10:
        print(f"  Runtime kisa ({runtime_min:.1f} dk) - buyuk batarya veya")
        print(f"  dusuk guc rejimi gerekli.")
    else:
        print(f"  Runtime yeterli ({runtime_min:.1f} dk).")

    if T_steady_default > COOLING_TARGET_T:
        print(f"  Pasif sogutma yetersiz (T = {T_steady_default:.1f} C > "
              f"{COOLING_TARGET_T} C)")
        print(f"  Aktif sogutma gerekli (h*A > {hA_needed:.2f} W/K)")
    else:
        print(f"  Pasif sogutma yeterli (T = {T_steady_default:.1f} C)")
    print()

    # --- Summary ---
    print("=" * 84)
    print("SONUC")
    print("=" * 84)
    print(f"  Power budget:  {P_tot:.2f} W")
    print(f"  Runtime:       {runtime_min:.1f} min")
    print(f"  Thermal load:  {q_density:.1f} W/m^2")
    print(f"  Cooling need:  {hA_needed:.2f} W/K")
    print()
    print("Not: Bu basit lumped-parameter model. Gercek thermal")
    print("simulasyon (FEM) sonraki asamalarda.")
    print("=" * 84)

    # --- Figure ---
    fig_path = Path("zt008_2_figure.png")
    fig_ok = save_thermal_figure(P_tot, str(fig_path))
    print(f"Figure saved: {fig_path} ({'OK' if fig_ok else 'FAIL'})")
    print()

    out = Path("zt008_2_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-008.2",
        "mode": a.mode,
        "power_breakdown_W": {k: float(v) for k, v in b.items()},
        "total_power_W": float(P_tot),
        "runtime_min": float(runtime_min),
        "thermal_steady_default_C": float(T_steady_default),
        "required_hA_W_per_K": float(hA_needed),
        "heat_flux_W_per_m2": float(q_density),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: lumped-parameter power + thermal only.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-008.2 tests."""
import numpy as np
from zt005.run_zt008_2 import (
    coil_resistance, joule_heating, power_breakdown, total_power,
    battery_runtime, thermal_steady_state, required_cooling_coeff,
    heat_flux_density, AMBIENT_T, COOLING_TARGET_T,
)


def test_coil_resistance_reasonable():
    R = coil_resistance()
    assert 0.01 < R < 1.0


def test_joule_heating_positive():
    P = joule_heating()
    assert P > 0


def test_power_breakdown_complete():
    b = power_breakdown()
    for key in ["coil_joule", "driver_loss", "control", "safety"]:
        assert key in b
        assert b[key] >= 0


def test_total_power_matches_sum():
    b = power_breakdown()
    assert abs(total_power() - sum(b.values())) < 1e-9


def test_battery_runtime_is_positive():
    r = battery_runtime(capacity_wh=100, power_w=50)
    assert r > 0


def test_battery_runtime_scales():
    r1 = battery_runtime(capacity_wh=100, power_w=50)
    r2 = battery_runtime(capacity_wh=100, power_w=100)
    assert abs(r1 / r2 - 2.0) < 1e-9


def test_thermal_steady_state_above_ambient():
    T = thermal_steady_state(10.0, cooling_coeff=1.0)
    assert T > AMBIENT_T


def test_required_cooling_for_target():
    hA = required_cooling_coeff(10.0, COOLING_TARGET_T)
    assert hA > 0


def test_heat_flux_positive():
    q = heat_flux_density(10.0)
    assert q > 0
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt008_2.py -v")
    print("  python -m zt005.run_zt008_2 --mode quick")