
"""ZT-008.4 - Sensor + telemetry design for backpack.

Covers: sensor selection, sampling rates, telemetry budget,
time synchronization strategy, and data logging capacity.
"""
import argparse, json
from pathlib import Path
import numpy as np


def sensor_list():
    """
    Sensor catalogue with type, range, sample rate, accuracy.
    """
    return [
        {
            "name": "3-axis B-field (Hall)",
            "quantity": "B",
            "range": "0 - 10 mT",
            "sample_rate_Hz": 10000,
            "accuracy_pct": 1.0,
            "interface": "SPI",
            "quantity_bytes": 6,
            "purpose": "Primary field measurement",
        },
        {
            "name": "3-axis E-field (capacitive)",
            "quantity": "E",
            "range": "0 - 100 V/m",
            "sample_rate_Hz": 10000,
            "accuracy_pct": 5.0,
            "interface": "analog (ADC)",
            "quantity_bytes": 6,
            "purpose": "Field consistency check",
        },
        {
            "name": "Coil current (Hall, DC-AC)",
            "quantity": "I",
            "range": "0 - 200 A",
            "sample_rate_Hz": 100000,
            "accuracy_pct": 0.5,
            "interface": "analog (ADC)",
            "quantity_bytes": 4,
            "purpose": "Safety + field scaling",
        },
        {
            "name": "Coil temperature (NTC)",
            "quantity": "T",
            "range": "-40 to 150 C",
            "sample_rate_Hz": 10,
            "accuracy_pct": 1.0,
            "interface": "analog (ADC)",
            "quantity_bytes": 4,
            "purpose": "Thermal safety",
        },
        {
            "name": "Battery voltage / current",
            "quantity": "V_bat, I_bat",
            "range": "0 - 60 V, -50 to +50 A",
            "sample_rate_Hz": 100,
            "accuracy_pct": 0.5,
            "interface": "I2C (BMS)",
            "quantity_bytes": 8,
            "purpose": "Power management",
        },
        {
            "name": "IMU (9-axis)",
            "quantity": "a, omega, mag",
            "range": "±16 g, ±2000 dps",
            "sample_rate_Hz": 1000,
            "accuracy_pct": 0.1,
            "interface": "SPI",
            "quantity_bytes": 18,
            "purpose": "Motion + safety",
        },
        {
            "name": "Reference clock (TCXO + GNSS)",
            "quantity": "t",
            "range": "atomic-disciplined",
            "sample_rate_Hz": 1,
            "accuracy_pct": 1e-10,   # 100 ps / s
            "interface": "UART + PPS",
            "quantity_bytes": 8,
            "purpose": "Time synchronization (critical)",
        },
        {
            "name": "Ambient environment (T, P, RH)",
            "quantity": "T_amb, P, RH",
            "range": "-40 to 85 C, 30-110 kPa, 0-100%",
            "sample_rate_Hz": 1,
            "accuracy_pct": 2.0,
            "interface": "I2C",
            "quantity_bytes": 6,
            "purpose": "Environmental baseline",
        },
    ]


def telemetry_budget():
    """
    Compute total data rate from all sensors (before compression).
    """
    total_bytes_per_sec = 0
    for s in sensor_list():
        # bytes per sample * sample rate
        rate = s["sample_rate_Hz"] * s["quantity_bytes"]
        total_bytes_per_sec += rate
    # Assume 30% compression from logging (delta encoding, etc.)
    return {
        "raw_bytes_per_s": total_bytes_per_sec,
        "compressed_bytes_per_s": int(total_bytes_per_sec * 0.7),
        "raw_MB_per_hour": total_bytes_per_sec * 3600 / 1e6,
    }


def channel_bandwidth():
    """BLE, WiFi, USB bandwidths (realistic effective)."""
    return {
        "BLE 5.0 (effective)": 100e3,       # 100 kB/s
        "WiFi 802.11n (effective)": 10e6,   # 10 MB/s
        "USB 2.0 (effective)": 30e6,        # 30 MB/s
        "USB 3.0 (effective)": 300e6,       # 300 MB/s
    }


def which_link_suffices(raw_rate_Bps):
    """Which telemetry link can carry the raw rate?"""
    result = {}
    for name, bw in channel_bandwidth().items():
        result[name] = bw >= raw_rate_Bps
    return result


def logging_capacity(runtime_hours, bytes_per_sec):
    """Storage required."""
    total_bytes = runtime_hours * 3600 * bytes_per_sec
    return total_bytes / 1e9  # GB


def time_sync_analysis():
    """
    Time synchronization strategy for causal validation.
    """
    return {
        "local_clock": {
            "type": "TCXO",
            "accuracy_ppm": 1.0,
            "drift_per_hour_us": 3600.0,
            "note": "Drifts fast - not sufficient alone",
        },
        "gnss_disciplined": {
            "type": "GNSS-disciplined TCXO",
            "accuracy_ppb": 10.0,
            "drift_per_hour_ns": 36.0,
            "note": "Requires GNSS lock - outdoor use",
        },
        "atomic_portable": {
            "type": "Chip-scale atomic clock (CSAC)",
            "accuracy_ppb": 0.01,
            "drift_per_hour_ns": 0.036,
            "note": "Best for indoor / mobile use",
        },
        "external_reference": {
            "type": "NTP/PTP over WiFi to lab clock",
            "accuracy_us": 1.0,
            "note": "Fallback / reference only",
        },
    }


def save_sensor_layout_figure(out_path):
    """Bar chart of sensor sample rates."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        sensors = sensor_list()
        names = [s["name"][:25] for s in sensors]
        rates = [s["sample_rate_Hz"] for s in sensors]

        fig, ax = plt.subplots(figsize=(11, 6))
        ax.barh(names, rates, color="#2ca02c")
        ax.set_xscale("log")
        ax.set_xlabel("Sample rate (Hz, log scale)")
        ax.set_title("Sensor sample rates")
        ax.grid(axis="x", alpha=0.3)
        plt.tight_layout()
        fig.savefig(out_path, dpi=120)
        plt.close(fig)
        return True
    except Exception as e:
        print(f"Figure save failed: {e}")
        return False


def save_telemetry_figure(out_path):
    """Bandwidth comparison."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        budget = telemetry_budget()
        links = channel_bandwidth()

        fig, ax = plt.subplots(figsize=(10, 6))
        names = list(links.keys()) + ["Sensor data (raw)", "Sensor data (comp.)"]
        vals = list(links.values()) + [
            budget["raw_bytes_per_s"],
            budget["compressed_bytes_per_s"],
        ]
        colors = ["#1f77b4"] * len(links) + ["#d62728", "#ff7f0e"]
        ax.bar(names, vals, color=colors)
        ax.set_yscale("log")
        ax.set_ylabel("Rate (bytes/s, log)")
        ax.set_title("Telemetry bandwidth vs sensor data rate")
        ax.grid(axis="y", alpha=0.3)
        plt.xticks(rotation=20, ha="right")
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
    ap.add_argument("--runtime-hours", type=float, default=2.0)
    a = ap.parse_args()

    print("=" * 100)
    print("ZAMANDA YOLCULUK - ZT-008.4")
    print("SENSOR + TELEMETRY DESIGN")
    print("=" * 100)
    print()

    # --- Sensor list ---
    print("=" * 100)
    print("SENSOR CATALOGUE")
    print("=" * 100)
    print(f"{'Sensor':<28s} {'Quantity':<12s} {'Range':<20s} "
          f"{'Rate (Hz)':>10s} {'Bytes':>7s}")
    print("-" * 100)
    for s in sensor_list():
        print(f"{s['name']:<28s} {s['quantity']:<12s} "
              f"{s['range']:<20s} {s['sample_rate_Hz']:>10.0f} "
              f"{s['quantity_bytes']:>7d}")
    print()

    # --- Telemetry budget ---
    print("=" * 100)
    print("TELEMETRY BUDGET")
    print("=" * 100)
    budget = telemetry_budget()
    print(f"  Raw data rate:        {budget['raw_bytes_per_s']:>15,d} bytes/s")
    print(f"  Compressed (30%):     {budget['compressed_bytes_per_s']:>15,d} bytes/s")
    print(f"  Raw per hour:         {budget['raw_MB_per_hour']:>15.2f} MB/h")
    print()

    print("  Link bandwidth:")
    for name, bw in channel_bandwidth().items():
        ok = "OK " if bw >= budget["raw_bytes_per_s"] else "NO "
        print(f"    [{ok}] {name:<30s} {bw:>15,.0f} bytes/s")
    print()

    # --- Logging capacity ---
    print("=" * 100)
    print("LOGGING CAPACITY")
    print("=" * 100)
    for hours in [0.5, 2.0, 8.0, 24.0]:
        cap = logging_capacity(hours, budget["compressed_bytes_per_s"])
        print(f"  {hours:>5.1f} h runtime -> {cap:.3f} GB compressed")
    print(f"  Recommended storage: 256 GB SSD (margins for multi-session)")
    print()

    # --- Time sync ---
    print("=" * 100)
    print("TIME SYNCHRONIZATION STRATEGY (critical for causal validation)")
    print("=" * 100)
    ts = time_sync_analysis()
    for name, data in ts.items():
        print(f"  {name}:")
        for k, v in data.items():
            print(f"    {k:<25s} {v}")
    print()

    print("  KARAR: Chip-scale atomic clock (CSAC) primary + GNSS/PTP fallback.")
    print("  Drift ~ 0.036 ns/hour -> 1 s in ~30,000 years.")
    print()

    # --- Figures ---
    fig1 = Path("zt008_4_sensors.png")
    fig2 = Path("zt008_4_telemetry.png")
    ok1 = save_sensor_layout_figure(str(fig1))
    ok2 = save_telemetry_figure(str(fig2))
    print(f"Figure saved: {fig1} ({'OK' if ok1 else 'FAIL'})")
    print(f"Figure saved: {fig2} ({'OK' if ok2 else 'FAIL'})")
    print()

    # --- Summary ---
    print("=" * 100)
    print("SONUC")
    print("=" * 100)
    print(f"  Sensors:              {len(sensor_list())}")
    print(f"  Raw data rate:        {budget['raw_bytes_per_s']:,} bytes/s")
    print(f"  Recommended link:     WiFi (10 MB/s) + USB fallback")
    print(f"  Storage (8h):         "
          f"{logging_capacity(8.0, budget['compressed_bytes_per_s']):.2f} GB")
    print(f"  Primary time source:  CSAC (chip-scale atomic clock)")
    print()
    print("Not: Bu kavramsal. Gercek secim icin EMC test, sensor")
    print("kalibrasyon, ve EN 55011 uyumlulugu gerekir.")
    print("=" * 100)

    out = Path("zt008_4_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-008.4",
        "mode": a.mode,
        "sensors": sensor_list(),
        "telemetry": budget,
        "logging_capacity_GB_per_8h":
            logging_capacity(8.0, budget["compressed_bytes_per_s"]),
        "time_sync": time_sync_analysis(),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: conceptual sensor + telemetry only.")


if __name__ == "__main__":
    main()
