
"""ZT-008.3 - Safety architecture for backpack device.

Layered safety: firmware -> MCU -> hardware -> mechanical.
All safety-critical functions must survive single-point failures.
"""
import argparse, json
from pathlib import Path
import numpy as np


def fault_modes():
    """
    FMEA-style table: each fault mode with severity, detection,
    and mitigation.
    """
    return [
        {
            "id": "F01",
            "component": "Coil overcurrent",
            "cause": "Driver failure / short circuit",
            "severity": "high",
            "detection": "Hall current sensor > 150 A",
            "mitigation": "Hardware relay cutoff < 10 ms",
            "layer": "hardware",
        },
        {
            "id": "F02",
            "component": "Coil overtemperature",
            "cause": "Cooling failure / excessive duty cycle",
            "severity": "high",
            "detection": "Thermistor > 70 C",
            "mitigation": "MCU cutoff + thermal fuse (T > 85 C)",
            "layer": "hardware + thermal",
        },
        {
            "id": "F03",
            "component": "Battery overvoltage / overcurrent",
            "cause": "BMS failure / short",
            "severity": "high",
            "detection": "BMS monitoring",
            "mitigation": "BMS cutoff + fuse",
            "layer": "hardware",
        },
        {
            "id": "F04",
            "component": "MCU firmware hang",
            "cause": "Software bug / EMI",
            "severity": "medium",
            "detection": "Hardware watchdog timer (500 ms)",
            "mitigation": "Reset MCU -> default SAFE state",
            "layer": "hardware",
        },
        {
            "id": "F05",
            "component": "Communication loss (BLE/WiFi)",
            "cause": "Radio interference / app crash",
            "severity": "low",
            "detection": "Timeout > 3 s (firmware) OR dead-man relay open",
            "mitigation": "Firmware auto-shutdown + hardware dead-man relay "
                          "(requires periodic heartbeat pulse)",
            "layer": "firmware + hardware",
        },
        {
            "id": "F06",
            "component": "Mechanical enclosure breach",
            "cause": "Drop / impact",
            "severity": "medium",
            "detection": "IMU spike > 5g (firmware) OR mechanical shock "
                         "switch trips at > 10g",
            "mitigation": "Firmware auto-shutdown + mechanical shock switch "
                          "hardwired to coil relay",
            "layer": "firmware + mechanical",
        },
        {
            "id": "F07",
            "component": "EMI interference with sensors",
            "cause": "Field resonance coupling",
            "severity": "medium",
            "detection": "Sensor consistency check",
            "mitigation": "Shielded sensor paths",
            "layer": "hardware",
        },
        {
            "id": "F08",
            "component": "Operator unauthorized start",
            "cause": "App compromise / accidental trigger",
            "severity": "high",
            "detection": "Two-factor: app + physical key",
            "mitigation": "Physical key required for ARM",
            "layer": "mechanical + firmware",
        },
    ]


def hardware_cutoff_chain():
    """
    Series of independent cutoff mechanisms.
    Any single mechanism can shut down the coil current.
    """
    return [
        {
            "name": "BMS cutoff",
            "trigger": "Battery OV/OC/OT",
            "response_ms": 20,
            "redundant": False,
        },
        {
            "name": "Current limit relay",
            "trigger": "I_coil > 150 A",
            "response_ms": 8,
            "redundant": True,
        },
        {
            "name": "Thermal fuse",
            "trigger": "T_coil > 85 C",
            "response_ms": 1000,
            "redundant": True,
        },
        {
            "name": "MCU interlock",
            "trigger": "Firmware state != RUN",
            "response_ms": 5,
            "redundant": True,
        },
        {
            "name": "Physical E-stop",
            "trigger": "Manual button",
            "response_ms": 2,
            "redundant": True,
        },
        {
            "name": "Watchdog reset",
            "trigger": "MCU hang > 500 ms",
            "response_ms": 500,
            "redundant": True,
        },
    ]


def state_machine():
    """
    Device state machine with transitions.
    """
    return {
        "OFF": {"ARM": "ARMED", "→": "SAFE"},
        "ARMED": {"START": "RUN", "ABORT": "OFF", "→": "SAFE"},
        "RUN":   {"STOP": "ARMED", "ABORT": "OFF",
                  "FAULT": "SAFE", "→": "SAFE"},
        "SAFE":  {"RESET": "OFF", "→": "OFF"},
    }


def check_safety_layers():
    """
    Verify each fault mode has at least one hardware-layer mitigation.
    Returns list of issues.
    """
    issues = []
    for fm in fault_modes():
        if "hardware" not in fm["layer"] and "mechanical" not in fm["layer"]:
            issues.append(f"{fm['id']}: no hardware/mechanical layer")
    return issues


def worst_case_response_time():
    """Worst-case response time across all cutoffs."""
    cutoffs = hardware_cutoff_chain()
    return max(c["response_ms"] for c in cutoffs)


def redundancy_analysis():
    """Count redundant vs single-point cutoffs."""
    cutoffs = hardware_cutoff_chain()
    redundant = sum(1 for c in cutoffs if c["redundant"])
    return {
        "total": len(cutoffs),
        "redundant": redundant,
        "single_point": len(cutoffs) - redundant,
    }


def save_safety_figure(out_path):
    """Plot cutoff response times."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        cutoffs = hardware_cutoff_chain()
        names = [c["name"] for c in cutoffs]
        times = [c["response_ms"] for c in cutoffs]

        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ["#2ca02c" if c["redundant"] else "#d62728"
                  for c in cutoffs]
        ax.barh(names, times, color=colors)
        ax.set_xscale("log")
        ax.set_xlabel("Response time (ms, log scale)")
        ax.set_title("Safety cutoff chain - response times")
        ax.axvline(10, color="k", linestyle="--", alpha=0.5,
                   label="Target: < 10 ms")
        ax.legend()
        ax.grid(axis="x", alpha=0.3)
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

    print("=" * 96)
    print("ZAMANDA YOLCULUK - ZT-008.3")
    print("SAFETY ARCHITECTURE")
    print("=" * 96)
    print()

    # --- FMEA table ---
    print("=" * 96)
    print("FMEA (Failure Mode and Effects Analysis)")
    print("=" * 96)
    print(f"{'ID':<5s} {'Component':<30s} {'Severity':<10s} "
          f"{'Detection':<30s} {'Layer':<18s}")
    print("-" * 96)
    for fm in fault_modes():
        print(f"{fm['id']:<5s} {fm['component']:<30s} "
              f"{fm['severity']:<10s} {fm['detection']:<30s} "
              f"{fm['layer']:<18s}")
    print()

    # --- Hardware cutoff chain ---
    print("=" * 96)
    print("HARDWARE CUTOFF CHAIN")
    print("=" * 96)
    print(f"{'Mechanism':<25s} {'Trigger':<30s} "
          f"{'Response (ms)':>15s} {'Redundant':>12s}")
    print("-" * 96)
    for c in hardware_cutoff_chain():
        red = "yes" if c["redundant"] else "NO"
        print(f"{c['name']:<25s} {c['trigger']:<30s} "
              f"{c['response_ms']:>15.0f} {red:>12s}")
    print()

    # --- State machine ---
    print("=" * 96)
    print("STATE MACHINE")
    print("=" * 96)
    for state, transitions in state_machine().items():
        print(f"  {state:<8s} -> {dict(transitions)}")
    print()

    # --- Analysis ---
    issues = check_safety_layers()
    wc_time = worst_case_response_time()
    red = redundancy_analysis()

    print("=" * 96)
    print("SAFETY ANALYSIS")
    print("=" * 96)
    print(f"  Fault modes:               {len(fault_modes())}")
    print(f"  Hardware/mechanical layers: "
          f"{len(fault_modes()) - len(issues)}/{len(fault_modes())}")
    print(f"  Issues:                    {len(issues)}")
    if issues:
        for i in issues:
            print(f"    - {i}")
    print()
    print(f"  Cutoff mechanisms:         {red['total']}")
    print(f"  Redundant:                 {red['redundant']}")
    print(f"  Single-point:              {red['single_point']}")
    print(f"  Worst-case response:       {wc_time:.0f} ms")
    print()

    # --- Assessment ---
    print("=" * 96)
    print("DEGERLENDIRME")
    print("=" * 96)
    if red["single_point"] == 0:
        print("  Tum cutoff mekanizmalari redundant.")
    else:
        print(f"  {red['single_point']} single-point cutoff var - "
              f"redundancy gerekli.")

    if wc_time > 10:
        print(f"  Worst-case response {wc_time:.0f} ms - hedefin üzerinde")
        print(f"  (hedef: < 10 ms)")
    else:
        print(f"  Worst-case response {wc_time:.0f} ms - hedef içinde")

    print()
    print("  Kritik kural: Yazilim guvenlik zincirinin TEK parcasi degildir.")
    print("  Her guvenlik fonksiyonu donanim veya mekanik katmanda")
    print("  tekrarlanmalidir.")
    print()

    # --- Summary ---
    print("=" * 96)
    print("SONUC")
    print("=" * 96)
    print(f"  FMEA:              {len(fault_modes())} modes")
    print(f"  Cutoff chain:      {red['total']} mechanisms "
          f"({red['redundant']} redundant)")
    print(f"  Worst response:    {wc_time:.0f} ms")
    print(f"  State machine:     {len(state_machine())} states")
    print()
    print("Not: Bu kavramsal. Gercek tasarimda IEC 61508 (SIL) veya")
    print("ISO 13849 (PL) standartlarina gore sertifikasyon gerekir.")
    print("=" * 96)

    # --- Figure ---
    fig_path = Path("zt008_3_figure.png")
    fig_ok = save_safety_figure(str(fig_path))
    print(f"Figure saved: {fig_path} ({'OK' if fig_ok else 'FAIL'})")
    print()

    out = Path("zt008_3_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-008.3",
        "mode": a.mode,
        "fault_modes": fault_modes(),
        "cutoff_chain": hardware_cutoff_chain(),
        "state_machine": state_machine(),
        "worst_case_response_ms": wc_time,
        "redundancy": red,
        "issues": issues,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: conceptual safety architecture only.")


if __name__ == "__main__":
    main()
