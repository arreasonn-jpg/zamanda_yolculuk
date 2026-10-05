
"""ZT-009.2 - Monte Carlo uncertainty propagation on h_00."""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.run_zt009_1 import h00_model, I_BASE, R_BASE, F_BASE, \
    ETA_BASE, D_BASE, H00_BASE

# Parameter uncertainties (relative, 1-sigma)
UNCERTAINTIES = {
    "I":   0.05,   # 5% current stability
    "R":   0.02,   # 2% mechanical tolerance
    "f":   0.001,  # 0.1% frequency stability (crystal)
    "eta": 0.10,   # 10% efficiency uncertainty
    "D":   0.01,   # 1% region positioning
}

DEFAULTS = {
    "I": I_BASE, "R": R_BASE, "f": F_BASE,
    "eta": ETA_BASE, "D": D_BASE,
}


def monte_carlo_h00(n_samples=10000, seed=42):
    """
    Sample each parameter from normal distribution with its
    relative uncertainty and propagate through h00_model.
    """
    rng = np.random.default_rng(seed)
    samples = {}
    for key, rel_sigma in UNCERTAINTIES.items():
        p0 = DEFAULTS[key]
        sigma = rel_sigma * p0
        samples[key] = rng.normal(p0, sigma, n_samples)

    # Clamp to physically valid ranges
    samples["eta"] = np.clip(samples["eta"], 0.0, 1.0)
    samples["I"]   = np.abs(samples["I"])
    samples["R"]   = np.abs(samples["R"])
    samples["D"]   = np.abs(samples["D"])

    h = h00_model(**samples)
    return h, samples


def one_at_a_time_uncertainty(n_samples=10000, seed=42):
    """
    Variance contribution from each parameter individually.
    """
    h_full, _ = monte_carlo_h00(n_samples, seed)
    var_full = np.var(h_full)

    contributions = {}
    for key in UNCERTAINTIES:
        # Vary only this parameter
        rng = np.random.default_rng(seed + 1)
        samples = {k: np.full(n_samples, DEFAULTS[k]) for k in DEFAULTS}
        p0 = DEFAULTS[key]
        sigma = UNCERTAINTIES[key] * p0
        samples[key] = rng.normal(p0, sigma, n_samples)
        if key == "eta":
            samples[key] = np.clip(samples[key], 0.0, 1.0)
        h = h00_model(**samples)
        contributions[key] = float(np.var(h) / var_full)
    return contributions


def summary_stats(h):
    return {
        "mean": float(np.mean(h)),
        "std": float(np.std(h)),
        "median": float(np.median(h)),
        "p05": float(np.percentile(h, 5)),
        "p95": float(np.percentile(h, 95)),
        "rel_uncertainty": float(np.std(h) / np.mean(h)),
    }


def save_histogram_figure(h, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        # Use log10(h) for visualization
        logh = np.log10(h)
        ax.hist(logh, bins=60, color="#1f77b4", alpha=0.7,
                edgecolor="black")
        ax.axvline(np.log10(H00_BASE), color="red", linestyle="--",
                   label=f"Baseline log10 = {np.log10(H00_BASE):.2f}")
        ax.set_xlabel("log10(h_00)")
        ax.set_ylabel("Frequency")
        ax.set_title(f"Monte Carlo uncertainty (N = {len(h)})")
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
    ap.add_argument("--n", type=int, default=10000)
    a = ap.parse_args()

    if a.mode == "full":
        a.n = 100000

    print("=" * 84)
    print("ZAMANDA YOLCULUK - ZT-009.2")
    print("MONTE CARLO UNCERTAINTY PROPAGATION")
    print("=" * 84)
    print(f"N samples = {a.n}")
    print()

    # --- Input uncertainties ---
    print("=" * 84)
    print("PARAMETER UNCERTAINTIES (1-sigma)")
    print("=" * 84)
    for key, sigma in UNCERTAINTIES.items():
        print(f"  {key:<6s}  {DEFAULTS[key]:.4e}  ±{sigma*100:.2f}%")
    print()

    # --- Monte Carlo ---
    h_samples, _ = monte_carlo_h00(a.n)
    stats = summary_stats(h_samples)

    print("=" * 84)
    print("OUTPUT DISTRIBUTION OF h_00")
    print("=" * 84)
    print(f"  Mean:            {stats['mean']:.4e}")
    print(f"  Median:          {stats['median']:.4e}")
    print(f"  Std:             {stats['std']:.4e}")
    print(f"  Relative std:    {stats['rel_uncertainty']*100:.3f}%")
    print(f"  90% CI:          [{stats['p05']:.4e}, {stats['p95']:.4e}]")
    print()

    # --- Variance contributions ---
    print("=" * 84)
    print("VARIANCE CONTRIBUTIONS")
    print("=" * 84)
    contrib = one_at_a_time_uncertainty(a.n)
    for key, frac in sorted(contrib.items(), key=lambda x: -x[1]):
        print(f"  {key:<6s}  {frac*100:>6.2f}%")
    print()

    # --- Compare sensitivity vs uncertainty ---
    print("=" * 84)
    print("SENSITIVITY vs UNCERTAINTY CONTRIBUTION")
    print("=" * 84)
    print(f"{'Param':<8s} {'|S_i|':>8s} {'RelSigma':>10s} "
          f"{'S*sigma':>12s} {'VarFrac':>10s}")
    print("-" * 84)
    from zt005.run_zt009_1 import sensitivity_index
    for key in UNCERTAINTIES:
        S = abs(sensitivity_index(key, rel_delta=0.001))
        rel = UNCERTAINTIES[key]
        product = S * rel
        frac = contrib[key]
        print(f"{key:<8s} {S:>8.3f} {rel:>10.4f} "
              f"{product:>12.4e} {frac:>10.2%}")
    print()

    # --- Interpretation ---
    print("=" * 84)
    print("DEGERLENDIRME")
    print("=" * 84)
    dominant = max(contrib.items(), key=lambda x: x[1])
    print(f"  Baskin parametre: {dominant[0]}  "
          f"({dominant[1]*100:.1f}% variyans)")
    print()
    print("  Bu, hangi parametrenin olcum/hardware belirsizliginin")
    print("  h_00 tahminine en cok katkida bulundugunu gosterir.")
    print()
    print("  Ancak baseline belirsizlik (relative_std) kucuk oldugu")
    print("  icin, gap 10^40 hala ana problem.")
    print()

    # --- Summary ---
    print("=" * 84)
    print("SONUC")
    print("=" * 84)
    print(f"  N:               {a.n}")
    print(f"  Mean h_00:       {stats['mean']:.4e}")
    print(f"  Rel uncertainty: {stats['rel_uncertainty']*100:.3f}%")
    print(f"  Dominant:        {dominant[0]} ({dominant[1]*100:.1f}%)")
    print()

    # --- Figure ---
    fig_path = Path("zt009_2_histogram.png")
    fig_ok = save_histogram_figure(h_samples, str(fig_path))
    print(f"Figure saved: {fig_path} ({'OK' if fig_ok else 'FAIL'})")
    print()

    out = Path("zt009_2_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-009.2",
        "mode": a.mode,
        "n_samples": a.n,
        "input_uncertainties": UNCERTAINTIES,
        "output_stats": stats,
        "variance_contributions": contrib,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: uncertainty propagation; no new physics.")


if __name__ == "__main__":
    main()
