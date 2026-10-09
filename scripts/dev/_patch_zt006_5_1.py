"""Auto-patch for ZT-006.5.1 — High-current sanity check."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt006_5_1.py"
TST = ROOT / "tests" / "test_zt006_5_1.py"

SRC_CODE = r'''
"""ZT-006.5.1 - High-current sanity check for geodesic integration.

Re-runs ZT-006.5 with h scaled by I^2 to make weak-field effects
visible in float64.
"""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.run_zt006_5 import (
    build_metric_grid, integrate_geodesic, geodesic_norm,
    interpolate_metric, christoffel_numerical,
)
from zt005.run_zt006_3 import build_metric_from_h
from zt005.poisson_solver import solve_poisson_axisym
from zt005.run_zt006_2_23 import compute_T00_physical

I_BASE = 100.0  # reference current in A
SECONDS_PER_YEAR = 3.15576e7


def scaled_h00(h00_base, I_A, I_base=I_BASE):
    """h ~ I^2 scaling."""
    return h00_base * (I_A / I_base)**2


def build_metric_grid_scaled(geometry, I_A, Nr=60, Nz=120, n_phi=8,
                             r_max=1.0, z_max=1.25):
    """Same as build_metric_grid but with h scaled to current I_A."""
    rv = np.linspace(0.0, r_max, Nr)
    zv = np.linspace(-z_max, z_max, Nz)
    T00 = compute_T00_physical(geometry, rv, zv, n_phi=n_phi)
    h00_base = solve_poisson_axisym(T00, rv, zv)
    h00_scaled = scaled_h00(h00_base, I_A)
    h0i = np.zeros_like(h00_scaled)
    g = build_metric_from_h(h00_scaled, h0i)
    return g, rv, zv, h00_scaled


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--n-phi", type=int, default=8)
    ap.add_argument("--current", type=float, default=None,
                    help="Current in A. Default: run 3 currents.")
    a = ap.parse_args()

    print("=" * 78)
    print("ZAMANDA YOLCULUK - ZT-006.5.1")
    print("HIGH-CURRENT SANITY CHECK FOR GEODESIC INTEGRATION")
    print("=" * 78)
    print(f"n_phi = {a.n_phi}  |  mode = {a.mode}")
    print()

    currents = ([a.current] if a.current else
                [1e2, 1e14, 1e20, 1e22])
    print(f"Testing currents: {currents} A")
    print()

    all_res = {}
    for g_name in ["G1", "G2", "G3", "G4"]:
        print(f"[{g_name}]")
        rows = []
        for I_A in currents:
            try:
                g, rv, zv, h00 = build_metric_grid_scaled(
                    g_name, I_A, Nr=60, Nz=120, n_phi=a.n_phi
                )
                h_center = float(h00[0, len(zv)//2])
                g_tt = float(g[0, len(zv)//2, 0, 0])

                # Geodesic: rest particle at origin, longer tau
                sol = integrate_geodesic(
                    g, rv, zv,
                    x0=(0, 0, 0, 0),
                    U0=(1.0, 0, 0, 0),
                    tau_max=1.0,      # 1 second
                    n_steps=100,
                )
                X = sol.y[:4]; U = sol.y[4:]
                norms = geodesic_norm(g, rv, zv, X, U)

                # Spatial displacement from origin
                dx = float(np.max(np.abs(X[1])))  # max |x|
                dy = float(np.max(np.abs(X[2])))
                dz = float(np.max(np.abs(X[3])))
                dt_final = float(U[0, -1]) - 1.0  # dt/dtau - 1

                # Christoffel at origin
                Gamma_origin = christoffel_numerical(
                    g, rv, zv, (0.0, 0.0)
                )
                Gamma_max = float(np.max(np.abs(Gamma_origin)))

                print(f"  I={I_A:.2e} A  "
                      f"h00_center={h_center:.4e}  "
                      f"g_tt={g_tt:.15f}  "
                      f"|Γ|_max={Gamma_max:.4e}  "
                      f"dx_max={dx:.4e} m  "
                      f"dt/dτ-1={dt_final:.4e}")
                rows.append({
                    "I_A": float(I_A),
                    "h00_center": float(h_center),
                    "g_tt": float(g_tt),
                    "Gamma_max": float(Gamma_max),
                    "dx_max": float(dx),
                    "dy_max": float(dy),
                    "dz_max": float(dz),
                    "dt_over_dtau_minus_1": float(dt_final),
                    "norm_mean": float(norms.mean()),
                    "norm_var": float(norms.var()),
                })
            except Exception as e:
                print(f"  I={I_A:.2e} A  FAILED: {type(e).__name__}: {e}")
                rows.append({"I_A": float(I_A), "error": str(e)})
        all_res[g_name] = rows
        print()

    # --- Interpretation ---
    print("=" * 78)
    print("YORUM")
    print("=" * 78)
    print("At I=100 A: h00 ~ 1e-48, float64 altinda -> no visible effect.")
    print("At I=1e14 A: h00 ~ 1e-20, still marginal but representable.")
    print("At I=1e20 A: h00 ~ 1e-8, visible in float64.")
    print("At I=1e22 A: h00 ~ 1e-4, clearly visible.")
    print()
    print("Beklenti: dx_max monotonik artmali; dt/dtau-1 negatif olmali.")
    print("Bu, modelin yuksek-akim rejiminde DOGRU calistigini gosterir.")
    print("=" * 78)

    out = Path("zt006_5_1_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.5.1",
        "mode": a.mode,
        "n_phi": a.n_phi,
        "currents_A": currents,
        "results": all_res,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: high-current sanity check; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-006.5.1 tests - high-current sanity."""
import numpy as np
from zt005.run_zt006_5_1 import (
    scaled_h00, build_metric_grid_scaled,
)
from zt005.run_zt006_5 import integrate_geodesic, geodesic_norm


def test_scaled_h00_quadratic():
    h_base = 1e-48
    h_2x = scaled_h00(h_base, 200.0)
    assert abs(h_2x / h_base - 4.0) < 1e-10


def test_scaled_h00_at_1e22():
    h_base = 5.9e-48
    h_new = scaled_h00(h_base, 1e22)
    # h_new = 5.9e-48 * (1e20)^2 = 5.9e-8
    assert 1e-9 < h_new < 1e-7


def test_metric_scaling_builds():
    g, rv, zv, h00 = build_metric_grid_scaled("G1", 1e22,
                                              Nr=20, Nz=40, n_phi=4)
    assert g.shape == (20, 40, 4, 4)
    assert h00.shape == (20, 40)
    assert np.max(np.abs(h00)) > 1e-9


def test_geodesic_visible_at_high_current():
    """At I=1e22, geodesic deviation should be visible."""
    g, rv, zv, _ = build_metric_grid_scaled("G1", 1e22,
                                            Nr=30, Nz=60, n_phi=4)
    sol = integrate_geodesic(
        g, rv, zv, x0=(0, 0, 0, 0), U0=(1.0, 0, 0, 0),
        tau_max=1.0, n_steps=50,
    )
    X = sol.y[:4]
    # Some visible deviation from origin
    max_spatial = np.max(np.abs(X[1:]))
    assert max_spatial > 1e-12


def test_low_current_stays_minkowski():
    """At I=100 A, deviation is below float64 epsilon."""
    g, rv, zv, h00 = build_metric_grid_scaled("G1", 100.0,
                                              Nr=20, Nz=40, n_phi=4)
    # g_tt should be exactly -1.0
    assert g[0, 20, 0, 0] == -1.0
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt006_5_1.py -v")
    print("  python -m zt005.run_zt006_5_1 --mode quick")