
"""ZT-008.6 - Prototype test protocol (safety-gated hierarchy)."""
import argparse, json
from pathlib import Path
import numpy as np


def test_hierarchy():
    """
    6-level test hierarchy. Each level requires previous level to pass.
    Human test is the last gate - requires all previous levels + IRB.
    """
    return [
        {
            "level": "T0",
            "name": "Simulation",
            "description": "Numerical verification (ZT-006.x pipeline)",
            "risk": "none",
            "prerequisites": [],
            "duration_h": 1,
            "human_involved": False,
            "gate_to_next": "All pipelines pass with rel_change < 1%",
        },
        {
            "level": "T1",
            "name": "Digital timestamp",
            "description": "Digital electronics + reference clock only. "
                           "Verify CSAC timing chain.",
            "risk": "none",
            "prerequisites": ["T0"],
            "duration_h": 4,
            "human_involved": False,
            "gate_to_next": "Clock drift < 1 ns/hour over 4 hours",
        },
        {
            "level": "T2",
            "name": "Shielded bench, low power",
            "description": "Coil + driver, current < 1 A, no field above "
                           "ambient. Validate control loop.",
            "risk": "low",
            "prerequisites": ["T1"],
            "duration_h": 8,
            "human_involved": False,
            "gate_to_next": "Control loop stable, no fault flags",
        },
        {
            "level": "T3",
            "name": "Bench with nominal power",
            "description": "Coil at 100 A, 1 MHz. Full sensor suite active. "
                           "Measure self-interference, thermal.",
            "risk": "medium",
            "prerequisites": ["T2"],
            "duration_h": 24,
            "human_involved": False,
            "gate_to_next": "Thermal < 40 C steady-state, EMI within limits",
        },
        {
            "level": "T4",
            "name": "Inorganic test object",
            "description": "Place inert mass in active volume. Verify no "
                           "anomalous effects (baseline for null hypothesis).",
            "risk": "medium",
            "prerequisites": ["T3"],
            "duration_h": 48,
            "human_involved": False,
            "gate_to_next": "Null result (no anomaly at >3 sigma)",
        },
        {
            "level": "T5",
            "name": "Biological model (non-human)",
            "description": "Cell cultures or tissue samples. Verify "
                           "safety envelope for biological matter.",
            "risk": "high",
            "prerequisites": ["T4"],
            "duration_h": 168,
            "human_involved": False,
            "gate_to_next": "Ethics committee + IRB approval",
        },
        {
            "level": "T6",
            "name": "Human test (LAST)",
            "description": "Only after all levels pass, plus independent "
                           "replication, plus legal/ethical review.",
            "risk": "extreme",
            "prerequisites": ["T5", "external_replication",
                              "ethics_committee", "irb",
                              "insurance", "regulatory"],
            "duration_h": 24,
            "human_involved": True,
            "gate_to_next": "NONE - end of protocol",
        },
    ]


def pretest_checklist(level):
    """Pre-test checklist for a given level."""
    base = [
        "Enclosure closed and labeled",
        "Emergency stop functional",
        "Watchdog timer verified",
        "Cooling path clear",
        "Reference clock locked",
    ]
    if level in ["T0", "T1"]:
        return base
    if level == "T2":
        return base + [
            "Current limit < 1 A verified",
            "RF shield in place",
            "Personnel clearance 1 m",
        ]
    if level == "T3":
        return base + [
            "Coil short-circuit test passed",
            "Thermal fuse verified",
            "Full sensor calibration current",
            "Personnel clearance 3 m",
        ]
    if level == "T4":
        return base + [
            "Test object mass verified",
            "Object position logged",
            "Independent observer present",
        ]
    if level == "T5":
        return base + [
            "Bio-safety cabinet ready",
            "Sample handling protocol reviewed",
            "Ethics committee acknowledgment on file",
        ]
    if level == "T6":
        return base + [
            "Medical personnel on standby",
            "Insurance activated",
            "Legal sign-off on file",
            "Independent safety officer present",
            "External replication confirmed",
        ]
    return base


def measurement_protocol(level):
    """Required measurements per level."""
    if level in ["T0", "T1"]:
        return {
            "duration_s": 3600,
            "sample_rate_hz": 1,
            "references": ["CSAC", "GNSS"],
            "metrics": ["clock_drift", "consistency"],
        }
    if level == "T2":
        return {
            "duration_s": 600,
            "sample_rate_hz": 10,
            "references": ["CSAC"],
            "metrics": ["current_stability", "temperature"],
        }
    if level == "T3":
        return {
            "duration_s": 3600,
            "sample_rate_hz": 10000,
            "references": ["CSAC", "GNSS"],
            "metrics": ["coil_current", "B_field", "E_field",
                        "temperature", "power"],
        }
    if level == "T4":
        return {
            "duration_s": 7200,
            "sample_rate_hz": 10000,
            "references": ["CSAC", "GNSS", "external_video"],
            "metrics": ["all_sensors", "object_position",
                        "electromagnetic_spectrum"],
        }
    return {
        "duration_s": 86400,
        "sample_rate_hz": 10000,
        "references": ["CSAC", "GNSS", "external_video",
                       "external_spectrum"],
        "metrics": ["all_sensors", "health_indicators"],
    }


def data_validation_criteria():
    """Run is valid if all criteria met."""
    return [
        {"criterion": "CSAC lock maintained",
         "test": "no unlock events"},
        {"criterion": "No safety faults",
         "test": "fault_count == 0"},
        {"criterion": "Temperature < 40 C",
         "test": "max_T < 40"},
        {"criterion": "Current within ±5%",
         "test": "|I - I_target| / I_target < 0.05"},
        {"criterion": "Sensor data continuous",
         "test": "no gaps > 100 ms"},
        {"criterion": "Data checksum valid",
         "test": "checksum == computed"},
    ]


def failure_response():
    """What to do when a test fails."""
    return [
        {"failure": "Clock unlock",
         "action": "Pause run, re-lock, restart from T1"},
        {"failure": "Safety fault",
         "action": "Full shutdown, root-cause analysis, restart from T2"},
        {"failure": "Thermal limit",
         "action": "Reduce duty cycle, re-verify T3 thermal model"},
        {"failure": "Anomalous measurement",
         "action": "Independent replication attempt; do NOT escalate"},
        {"failure": "Data corruption",
         "action": "Redo run; no partial credit"},
    ]


def human_gate_requirements():
    """
    6-layer gate for T6 (human test). All must pass.
    """
    return [
        {"layer": 1, "name": "Physics",
         "requirement": "T0-T5 all pass + external replication"},
        {"layer": 2, "name": "Engineering",
         "requirement": "Hardware safety certified (SIL 2+)"},
        {"layer": 3, "name": "Biology",
         "requirement": "Bio-model tests show no harm at 10x dose"},
        {"layer": 4, "name": "Ethics",
         "requirement": "Ethics committee + IRB approval"},
        {"layer": 5, "name": "Legal",
         "requirement": "Insurance + regulatory compliance"},
        {"layer": 6, "name": "Medical",
         "requirement": "Medical team + emergency response plan"},
    ]


def save_hierarchy_figure(out_path):
    """Plot test hierarchy as pyramid."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        levels = test_hierarchy()
        n = len(levels)

        fig, ax = plt.subplots(figsize=(10, 8))
        colors = ["#2ca02c", "#2ca02c", "#f5d76e", "#f5d76e",
                  "#f39c12", "#e67e22", "#c0392b"]

        for i, lvl in enumerate(levels):
            y = i
            x_center = 0.5
            width = 0.4 + 0.6 * (i / (n - 1))
            ax.add_patch(plt.Rectangle(
                (x_center - width/2, y - 0.4), width, 0.8,
                facecolor=colors[i], edgecolor="black", alpha=0.7
            ))
            label = f"{lvl['level']}: {lvl['name']}"
            if lvl["human_involved"]:
                label += " [HUMAN]"
            ax.text(x_center, y, label, ha="center", va="center",
                    fontsize=10, fontweight="bold" if lvl["human_involved"]
                    else "normal")

        ax.set_xlim(0, 1)
        ax.set_ylim(-1, n)
        ax.axis("off")
        ax.set_title("Test hierarchy (safety-gated pyramid)",
                     fontsize=13)
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

    print("=" * 100)
    print("ZAMANDA YOLCULUK - ZT-008.6")
    print("PROTOTYPE TEST PROTOCOL")
    print("=" * 100)
    print()

    hierarchy = test_hierarchy()

    print("=" * 100)
    print("TEST HIERARCHY")
    print("=" * 100)
    print(f"{'Level':<6s} {'Name':<28s} {'Risk':<10s} "
          f"{'Human':<7s} {'Duration':<12s} Prereq")
    print("-" * 100)
    for t in hierarchy:
        h = "yes" if t["human_involved"] else "no"
        prereq = ",".join(t["prerequisites"]) if t["prerequisites"] else "-"
        print(f"{t['level']:<6s} {t['name']:<28s} {t['risk']:<10s} "
              f"{h:<7s} {t['duration_h']:>5.0f} h      {prereq}")
    print()

    # --- Pre-test checklist per level ---
    print("=" * 100)
    print("PRE-TEST CHECKLIST SUMMARY")
    print("=" * 100)
    for t in hierarchy:
        items = pretest_checklist(t["level"])
        print(f"  {t['level']}  ({len(items)} items):")
        for it in items[:3]:
            print(f"    - {it}")
        if len(items) > 3:
            print(f"    ... (+{len(items) - 3} more)")
    print()

    # --- Measurement protocol ---
    print("=" * 100)
    print("MEASUREMENT PROTOCOL")
    print("=" * 100)
    print(f"{'Level':<6s} {'Duration':<12s} {'Rate (Hz)':<12s} "
          f"{'References':<30s} Metrics")
    print("-" * 100)
    for t in hierarchy:
        mp = measurement_protocol(t["level"])
        refs = ",".join(mp["references"])
        metrics = ",".join(mp["metrics"])[:40]
        print(f"{t['level']:<6s} {mp['duration_s']:>8.0f} s  "
              f"{mp['sample_rate_hz']:>8.0f}    {refs:<30s} {metrics}")
    print()

    # --- Data validation ---
    print("=" * 100)
    print("DATA VALIDATION CRITERIA")
    print("=" * 100)
    for c in data_validation_criteria():
        print(f"  [ ] {c['criterion']:<35s} ({c['test']})")
    print()

    # --- Failure response ---
    print("=" * 100)
    print("FAILURE RESPONSE MATRIX")
    print("=" * 100)
    for f in failure_response():
        print(f"  {f['failure']:<25s} -> {f['action']}")
    print()

    # --- Human gate ---
    print("=" * 100)
    print("HUMAN TEST GATE (T6) - 6 LAYERS")
    print("=" * 100)
    for g in human_gate_requirements():
        print(f"  Layer {g['layer']}: {g['name']:<15s} "
              f"- {g['requirement']}")
    print()

    # --- Summary ---
    print("=" * 100)
    print("SONUC")
    print("=" * 100)
    total_hours = sum(t["duration_h"] for t in hierarchy)
    human_levels = sum(1 for t in hierarchy if t["human_involved"])
    print(f"  Test levels:         {len(hierarchy)}")
    print(f"  Total planned hours: {total_hours}")
    print(f"  Human-involving:     {human_levels}")
    print(f"  Prerequisite chain:  T0 -> T1 -> ... -> T6")
    print()
    print("  KRITIK: Insan deneyi (T6) tek basina 'bir sonraki test'")
    print("  degildir. 6 katmanli gate + harici replikasyon +")
    print("  etik/hukuki/tibbi onay gerektirir.")
    print()

    fig_path = Path("zt008_6_hierarchy.png")
    fig_ok = save_hierarchy_figure(str(fig_path))
    print(f"Figure saved: {fig_path} ({'OK' if fig_ok else 'FAIL'})")
    print()

    out = Path("zt008_6_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-008.6",
        "mode": a.mode,
        "test_hierarchy": hierarchy,
        "data_validation": data_validation_criteria(),
        "failure_response": failure_response(),
        "human_gate": human_gate_requirements(),
        "total_planned_hours": total_hours,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: conceptual test protocol only.")


if __name__ == "__main__":
    main()
