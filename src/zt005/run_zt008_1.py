
"""ZT-008.1 - Backpack mechanical layout + component BOM."""
import argparse, json
from pathlib import Path
import numpy as np

# Hard envelope (ZT-002)
ENVELOPE_CM = (45.0, 32.0, 15.0)   # H x W x D
ENVELOPE_M3 = (0.45, 0.32, 0.15)


def component_list():
    """
    List of components with dimensions (cm), mass (kg), and
    placement (fractional position within envelope).
    """
    return [
        # name, L_cm, W_cm, H_cm, mass_kg, (x, y, z) center fracs
        {
            "name": "EM coil array (6 loops)",
            "dims_cm": (12.0, 12.0, 12.0),
            "mass_kg": 2.5,
            "position_frac": (0.10, 0.50, 0.50),
            "category": "field",
        },
        {
            "name": "HV driver + resonant tank",
            "dims_cm": (10.0, 8.0, 6.0),
            "mass_kg": 1.2,
            "position_frac": (0.35, 0.30, 0.50),
            "category": "power",
        },
        {
            "name": "Power module (battery + BMS)",
            "dims_cm": (12.0, 10.0, 5.0),
            "mass_kg": 2.0,
            "position_frac": (0.60, 0.50, 0.50),
            "category": "power",
        },
        {
            "name": "Control MCU (ESP32-S3)",
            "dims_cm": (6.0, 4.0, 1.5),
            "mass_kg": 0.05,
            "position_frac": (0.80, 0.40, 0.50),
            "category": "control",
        },
        {
            "name": "SBC (Raspberry Pi CM4)",
            "dims_cm": (8.0, 6.0, 2.0),
            "mass_kg": 0.15,
            "position_frac": (0.80, 0.70, 0.50),
            "category": "control",
        },
        {
            "name": "EM field sensors",
            "dims_cm": (4.0, 3.0, 1.0),
            "mass_kg": 0.05,
            "position_frac": (0.50, 0.15, 0.30),
            "category": "sensor",
        },
        {
            "name": "Thermal + IMU sensors",
            "dims_cm": (4.0, 3.0, 1.0),
            "mass_kg": 0.05,
            "position_frac": (0.50, 0.85, 0.30),
            "category": "sensor",
        },
        {
            "name": "Safety relay + e-stop",
            "dims_cm": (5.0, 4.0, 3.0),
            "mass_kg": 0.30,
            "position_frac": (0.90, 0.50, 0.50),
            "category": "safety",
        },
        {
            "name": "Cooling + thermal paths",
            "dims_cm": (30.0, 2.0, 1.0),
            "mass_kg": 0.40,
            "position_frac": (0.50, 0.50, 0.85),
            "category": "thermal",
        },
        {
            "name": "Structure + mounting frame",
            "dims_cm": (44.0, 31.0, 14.0),
            "mass_kg": 1.5,
            "position_frac": (0.50, 0.50, 0.50),
            "category": "structure",
        },
    ]


def check_fit(components, envelope_cm=ENVELOPE_CM):
    """
    Check each component fits within envelope.
    Return list of (component_name, ok, note).
    """
    results = []
    H, W, D = envelope_cm
    for c in components:
        L, w, h = c["dims_cm"]
        fits = (L <= H and w <= W and h <= D)
        results.append({
            "name": c["name"],
            "dims_cm": (L, w, h),
            "fits": fits,
        })
    return results


def total_mass(components):
    return sum(c["mass_kg"] for c in components)


def category_mass(components):
    cats = {}
    for c in components:
        cats[c["category"]] = cats.get(c["category"], 0.0) + c["mass_kg"]
    return cats


def print_layout(components):
    H, W, D = ENVELOPE_CM
    print(f"Envelope: {H} x {W} x {D} cm  (H x W x D)")
    print()
    print(f"{'Component':<30s} {'Dims (cm)':<18s} "
          f"{'Mass (kg)':>10s} {'Category':<12s} "
          f"{'Pos (frac)':<18s}")
    print("-" * 96)
    for c in components:
        L, w, h = c["dims_cm"]
        pos = c["position_frac"]
        print(f"{c['name']:<30s} {f'{L:.0f}x{w:.0f}x{h:.0f}':<18s} "
              f"{c['mass_kg']:>10.3f} {c['category']:<12s} "
              f"({pos[0]:.2f},{pos[1]:.2f},{pos[2]:.2f})")
    print()


def print_fit_check(components):
    results = check_fit(components)
    print("HARD CONSTRAINT CHECK (each component vs envelope):")
    print()
    all_ok = True
    for r in results:
        status = "OK " if r["fits"] else "FAIL"
        if not r["fits"]:
            all_ok = False
        L, w, h = r["dims_cm"]
        print(f"  [{status}] {r['name']:<30s} {L:.0f}x{w:.0f}x{h:.0f} cm")
    print()
    print(f"All components fit: {all_ok}")
    return all_ok


def save_layout_figure(components, out_path):
    """Simple 3D scatter of component centers."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa

        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection="3d")
        H, W, D = ENVELOPE_M3

        # Envelope wireframe corners
        corners = np.array([
            [0, 0, 0], [H, 0, 0], [H, W, 0], [0, W, 0],
            [0, 0, D], [H, 0, D], [H, W, D], [0, W, D],
        ])
        edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
                 (0,4),(1,5),(2,6),(3,7)]
        for a, b in edges:
            ax.plot(*zip(corners[a], corners[b]), "k-", alpha=0.3, linewidth=0.8)

        for c in components:
            xf, yf, zf = c["position_frac"]
            x = xf * H; y = yf * W; z = zf * D
            ax.scatter(x, y, z, s=c["mass_kg"] * 100 + 20,
                       label=c["name"][:20])
            ax.text(x, y, z, c["name"][:12], fontsize=7)

        ax.set_xlabel("Height (m)")
        ax.set_ylabel("Width (m)")
        ax.set_zlabel("Depth (m)")
        ax.set_title("Backpack component layout (45 x 32 x 15 cm)")
        ax.legend(loc="upper left", fontsize=7)
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
    print("ZAMANDA YOLCULUK - ZT-008.1")
    print("BACKPACK MECHANICAL LAYOUT")
    print("=" * 96)
    print()

    comps = component_list()
    print_layout(comps)
    ok = print_fit_check(comps)

    mass_total = total_mass(comps)
    cats = category_mass(comps)

    print("=" * 96)
    print("MASS BUDGET")
    print("=" * 96)
    for cat, m in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {cat:<15s} {m:.3f} kg")
    print("-" * 40)
    print(f"  {'TOTAL':<15s} {mass_total:.3f} kg")
    print()
    print(f"  Target: backpack should weigh < 8 kg")
    print(f"  Actual: {mass_total:.2f} kg")
    print(f"  Headroom: {8.0 - mass_total:.2f} kg")
    print()

    # --- Save figure ---
    fig_path = Path("zt008_1_figure.png")
    fig_ok = save_layout_figure(comps, str(fig_path))
    print(f"Figure saved: {fig_path} ({'OK' if fig_ok else 'FAIL'})")
    print()

    # --- Summary ---
    print("=" * 96)
    print("SONUC")
    print("=" * 96)
    print(f"  Components: {len(comps)}")
    print(f"  Envelope:   {ENVELOPE_CM[0]} x {ENVELOPE_CM[1]} x {ENVELOPE_CM[2]} cm")
    print(f"  Total mass: {mass_total:.2f} kg")
    print(f"  Fit check:  {'PASS' if ok else 'FAIL'}")
    print()
    print("Not: Bu layout kavramsal. Gercek parca secimi ve")
    print("detayli CAE (thermal, structural, EMI) sonraki asamalarda.")
    print("=" * 96)

    out = Path("zt008_1_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-008.1",
        "mode": a.mode,
        "envelope_cm": ENVELOPE_CM,
        "components": comps,
        "total_mass_kg": mass_total,
        "category_mass": cats,
        "fit_ok": ok,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: conceptual layout only.")


if __name__ == "__main__":
    main()
