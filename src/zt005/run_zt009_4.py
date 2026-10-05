
"""ZT-009.4 - Generate final paper figures."""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.run_zt007_b7 import mechanism_table
from zt005.run_zt009_1 import (
    h00_model, sensitivity_index, H00_BASE,
    I_BASE, R_BASE, F_BASE, ETA_BASE, D_BASE,
)
from zt005.run_zt009_2 import monte_carlo_h00, UNCERTAINTIES
from zt005.run_zt009_3 import (
    tier_scenarios, compute_h00_with_engineering, H00_TARGET,
)
from zt005.run_zt008_1 import component_list

FIGDIR = Path("paper_figures")


def fig1_unified_mechanisms(outdir):
    """Figure 1: 6 mechanisms unified bar chart."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        mechs = mechanism_table()
        names = [m["key"] for m in mechs]
        log_gaps = [np.log10(m["gap"]) for m in mechs]
        colors = ["#2ca02c" if m["physically_allowed"] else "#d62728"
                  for m in mechs]

        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.barh(names, log_gaps, color=colors)
        ax.set_xlabel("log10 (Required / Device)", fontsize=12)
        ax.set_title("Figure 1: Feasibility gap across 6 mechanisms",
                     fontsize=13)
        ax.set_xlim(0, max(log_gaps) + 5)
        ax.grid(axis="x", alpha=0.3)
        for bar, lg in zip(bars, log_gaps):
            ax.text(bar.get_width() + 0.5,
                    bar.get_y() + bar.get_height()/2,
                    f"10^{lg:.1f}", va="center", fontsize=10)
        plt.tight_layout()
        path = outdir / "fig1_unified_mechanisms.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path
    except Exception as e:
        print(f"Fig1 failed: {e}")
        return None


def fig2_em_scaling(outdir):
    """Figure 2: h_00 vs current (scaling study)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        currents = np.logspace(2, 22, 41)  # 100 A to 10^22 A
        h_values = np.array([h00_model(I=I) for I in currents])
        target = H00_TARGET

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.loglog(currents, h_values, "b-", linewidth=2,
                  label=r"$h_{00}(I)$ (this work)")
        ax.axhline(target, color="r", linestyle="--",
                   label=f"Target (1 s/yr)")
        ax.axvline(1e22, color="k", linestyle=":",
                   label="Required I ~ 10^22 A")

        # Mark reference currents
        for I_ref, name in [(1e5, "Lightning"),
                             (1e7, "Lab pulse"),
                             (1e18, "Astro"),
                             (1e24, "Magnetar")]:
            ax.axvline(I_ref, color="gray", alpha=0.3, linewidth=0.5)
            ax.text(I_ref, target * 1e-15, name, rotation=90,
                    fontsize=8, va="bottom")

        ax.set_xlabel("Current I (A)", fontsize=12)
        ax.set_ylabel(r"$h_{00}$", fontsize=12)
        ax.set_title("Figure 2: EM-induced metric perturbation scaling",
                     fontsize=13)
        ax.legend(loc="lower right")
        ax.grid(True, which="both", alpha=0.3)
        plt.tight_layout()
        path = outdir / "fig2_em_scaling.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path
    except Exception as e:
        print(f"Fig2 failed: {e}")
        return None


def fig3_best_case_tiers(outdir):
    """Figure 3: Best-case engineering tiers."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        tiers = tier_scenarios()
        names = [t["name"].split("(")[0].strip() for t in tiers]
        h_values = []
        for t in tiers:
            h = compute_h00_with_engineering(
                t["I_A"], t["R_m"], t["eta"],
                t["N_stages"], t["Q_factor"],
            )
            h_values.append(h)

        log_h = np.log10(h_values)
        log_target = np.log10(H00_TARGET)
        log_base = np.log10(H00_BASE)

        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ["#999999", "#f39c12", "#2ca02c"]
        bars = ax.barh(names, log_h, color=colors)

        ax.axvline(log_target, color="red", linestyle="--",
                   label=f"Target: 10^{log_target:.1f}")
        ax.axvline(log_base, color="blue", linestyle=":",
                   label=f"Baseline: 10^{log_base:.1f}")

        for bar, lg in zip(bars, log_h):
            ax.text(lg - 1, bar.get_y() + bar.get_height()/2,
                    f"10^{lg:.1f}", va="center", ha="right",
                    fontsize=11, color="white", fontweight="bold")

        ax.set_xlabel("log10(h_00)", fontsize=12)
        ax.set_title("Figure 3: Best-case engineering scenarios",
                     fontsize=13)
        ax.set_xlim(log_base - 10, log_target + 5)
        ax.legend(loc="lower right")
        ax.grid(axis="x", alpha=0.3)
        plt.tight_layout()
        path = outdir / "fig3_best_case_tiers.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path
    except Exception as e:
        print(f"Fig3 failed: {e}")
        return None


def fig4_monte_carlo(outdir):
    """Figure 4: Monte Carlo histogram."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        h, _ = monte_carlo_h00(50000)
        log_h = np.log10(h)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(log_h, bins=80, color="#1f77b4", alpha=0.75,
                edgecolor="black", density=True)
        ax.axvline(np.log10(H00_BASE), color="red", linestyle="--",
                   linewidth=2,
                   label=f"Mean log10 = {np.log10(H00_BASE):.2f}")
        ax.set_xlabel("log10(h_00)", fontsize=12)
        ax.set_ylabel("Probability density", fontsize=12)
        ax.set_title(f"Figure 4: Monte Carlo (N = 50000) - relative std = {np.std(h)/np.mean(h)*100:.1f}%",
                     fontsize=12)
        ax.legend()
        ax.grid(alpha=0.3)
        plt.tight_layout()
        path = outdir / "fig4_monte_carlo.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path
    except Exception as e:
        print(f"Fig4 failed: {e}")
        return None


def fig5_sensitivity_tornado(outdir):
    """Figure 5: Sensitivity tornado chart."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        params = ["I", "R", "f", "eta", "D"]
        labels = ["Current I", "Radius R", "Frequency f",
                  "Efficiency η", "Diameter D"]
        S_values = [sensitivity_index(p, rel_delta=0.01) for p in params]

        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ["#2ca02c" if s >= 0 else "#d62728" for s in S_values]
        y_pos = np.arange(len(params))
        ax.barh(y_pos, S_values, color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel(r"$S_i = (dh/h)/(dp/p)$", fontsize=12)
        ax.set_title("Figure 5: Sensitivity of h_00 to device parameters",
                     fontsize=13)
        ax.axvline(0, color="k", linewidth=0.5)
        for i, s in enumerate(S_values):
            ax.text(s + 0.05 * np.sign(s), i, f"{s:+.2f}",
                    va="center", ha="left" if s > 0 else "right",
                    fontsize=10)
        ax.grid(axis="x", alpha=0.3)
        plt.tight_layout()
        path = outdir / "fig5_sensitivity_tornado.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path
    except Exception as e:
        print(f"Fig5 failed: {e}")
        return None


def fig6_engineering_envelope(outdir):
    """Figure 6: Backpack component layout."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa

        comps = component_list()
        H, W, D = 0.45, 0.32, 0.15

        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection="3d")

        corners = np.array([
            [0, 0, 0], [H, 0, 0], [H, W, 0], [0, W, 0],
            [0, 0, D], [H, 0, D], [H, W, D], [0, W, D],
        ])
        edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
                 (0,4),(1,5),(2,6),(3,7)]
        for a, b in edges:
            ax.plot(*zip(corners[a], corners[b]), "k-", alpha=0.3,
                    linewidth=0.8)

        for c in comps:
            xf, yf, zf = c["position_frac"]
            x = xf * H; y = yf * W; z = zf * D
            ax.scatter(x, y, z, s=c["mass_kg"] * 120 + 20,
                       label=c["name"][:20])
            ax.text(x, y, z, c["name"][:10], fontsize=7)

        ax.set_xlabel("Height (m)", fontsize=10)
        ax.set_ylabel("Width (m)", fontsize=10)
        ax.set_zlabel("Depth (m)", fontsize=10)
        ax.set_title("Figure 6: Backpack component layout",
                     fontsize=13)
        plt.tight_layout()
        path = outdir / "fig6_backpack_layout.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path
    except Exception as e:
        print(f"Fig6 failed: {e}")
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    FIGDIR.mkdir(exist_ok=True)

    print("=" * 84)
    print("ZAMANDA YOLCULUK - ZT-009.4")
    print("FINAL PAPER FIGURES")
    print("=" * 84)
    print(f"Output directory: {FIGDIR}")
    print()

    figures = []
    for i, fn in enumerate([fig1_unified_mechanisms,
                             fig2_em_scaling,
                             fig3_best_case_tiers,
                             fig4_monte_carlo,
                             fig5_sensitivity_tornado,
                             fig6_engineering_envelope], 1):
        print(f"[{i}/6] Generating figure {i}...")
        path = fn(FIGDIR)
        if path:
            figures.append(str(path))
            print(f"      -> {path}")
        else:
            print(f"      FAILED")
    print()

    print("=" * 84)
    print("SONUC")
    print("=" * 84)
    print(f"  Figures generated: {len(figures)}/6")
    for f in figures:
        print(f"    {f}")
    print()

    # --- Summary manifest ---
    manifest = {
        "stage": "ZT-009.4",
        "mode": a.mode,
        "figures": figures,
        "figure_descriptions": {
            "fig1": "6 mechanisms unified bar chart",
            "fig2": "EM-induced metric perturbation scaling",
            "fig3": "Best-case engineering scenarios",
            "fig4": "Monte Carlo uncertainty histogram",
            "fig5": "Sensitivity tornado chart",
            "fig6": "Backpack component layout",
        },
    }
    out = Path("zt009_4_results.json")
    out.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Manifest written to {out}")
    print("SCIENTIFIC STATUS: paper figure generation complete.")


if __name__ == "__main__":
    main()
