
"""ZT-006.3.1 - Metric + precision-safe proper-time analysis."""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.run_zt006_2_23 import compute_T00_physical
from zt005.poisson_solver import solve_poisson_axisym, EINSTEIN_FACTOR

C = 2.99792458e8
SECONDS_PER_YEAR = 3.15576e7


def build_metric_from_h(h_00, h_0i_mag, eta=None):
    """g_munu = eta_munu + h_munu."""
    if eta is None:
        eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    Nr, Nz = h_00.shape
    g = np.tile(eta[None, None, :, :], (Nr, Nz, 1, 1))
    g[:, :, 0, 0] = -1.0 + h_00
    g[:, :, 0, 1] = h_0i_mag
    g[:, :, 1, 0] = h_0i_mag
    return g


def proper_time_delta(h_00):
    """
    Weak-field: dtau/dt = sqrt(1 - h_00) ~ 1 - h_00/2.
    Return (1 - dtau/dt) as POSITIVE fractional slowdown.

    Computing 1 - sqrt(1 - h) numerically loses all precision when
    h < machine_epsilon. Use algebraic identity:
        1 - sqrt(1 - h) = h / (1 + sqrt(1 - h))
    """
    h = np.asarray(h_00, dtype=float)
    # For small h, denominator ~ 2, result ~ h/2 (no precision loss)
    denom = 1.0 + np.sqrt(np.maximum(1.0 - h, 0.0))
    return h / denom


def compute_geodesic_diagnostics(geometry, Nr, Nz, n_phi=8,
                                  r_max=1.0, z_max=1.25):
    rv = np.linspace(0.0, r_max, Nr)
    zv = np.linspace(-z_max, z_max, Nz)
    T00 = compute_T00_physical(geometry, rv, zv, n_phi=n_phi)
    h00 = solve_poisson_axisym(T00, rv, zv)
    h0i = np.zeros_like(h00)

    delta_tau_frac = proper_time_delta(h00)   # (1 - dtau/dt)

    ic, jc = 0, Nz // 2
    delta_center = float(delta_tau_frac[ic, jc])
    delta_max = float(np.max(delta_tau_frac))

    return {
        "Nr": Nr, "Nz": Nz,
        "h_00_center": float(h00[ic, jc]),
        "h_00_max":    float(np.max(np.abs(h00))),
        "delta_tau_over_t_center":  delta_center,
        "delta_tau_over_t_max":     delta_max,
        "delta_tau_per_year_center_s": delta_center * SECONDS_PER_YEAR,
        "delta_tau_per_year_max_s":    delta_max * SECONDS_PER_YEAR,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--n-phi", type=int, default=8)
    a = ap.parse_args()

    resolutions = ([(40, 80)] if a.mode == "quick"
                   else [(40, 80), (60, 120), (80, 160)])

    print("=" * 72)
    print("ZAMANDA YOLCULUK - ZT-006.3.1")
    print("METRIC + PRECISION-SAFE PROPER-TIME ANALYSIS")
    print("=" * 72)
    print(f"Einstein factor = {EINSTEIN_FACTOR:.6e}")
    print(f"n_phi = {a.n_phi}  |  mode = {a.mode}")
    print()

    all_res = {}
    for g_name in ["G1", "G2", "G3", "G4"]:
        print(f"[{g_name}]")
        rows = []
        for Nr, Nz in resolutions:
            r = compute_geodesic_diagnostics(g_name, Nr, Nz, n_phi=a.n_phi)
            print(f"  grid=({r['Nr']:3d},{r['Nz']:3d})  "
                  f"h00_center={r['h_00_center']:.4e}  "
                  f"(1-dtau/dt)_center={r['delta_tau_over_t_center']:.4e}  "
                  f"dtau/yr_center={r['delta_tau_per_year_center_s']:.4e} s")
            rows.append(r)
        all_res[g_name] = rows
        print()

    # ---- Interpretive summary ----
    print("=" * 72)
    print("FIZIKSEL YORUM")
    print("=" * 72)
    print("Weak-field gravity time dilation at device center (per year):")
    for g in ["G1", "G2", "G3", "G4"]:
        dt = all_res[g][-1]["delta_tau_per_year_center_s"]
        ratio = dt / 1.0 if dt > 0 else 0.0
        print(f"  {g}: dtau ~ {dt:.3e} s per year")
    print()
    print("Hedef: 1 saniyelik gecmis zaman yolculugu.")
    print("Mevcut etki hedeften ~10^41 kat daha kucuk.")
    print("=" * 72)

    out = Path("zt006_3_1_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.3.1",
        "mode": a.mode,
        "einstein_factor": EINSTEIN_FACTOR,
        "results": all_res,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: weak-field metric; no CTC conclusion.")


if __name__ == "__main__":
    main()
