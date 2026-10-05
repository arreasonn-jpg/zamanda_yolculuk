"""ZT-006.2.21 runner: Poisson-solver convergence study."""
import argparse
import json
from pathlib import Path
import numpy as np

from zt005.poisson_solver import solve_poisson_axisym, EINSTEIN_FACTOR


def make_test_source(r_vals, z_vals, geometry="G1"):
    """
    Synthetic T00 source for numerical convergence tests.
    Real source builders (from ZT-006.2.x) will be plugged in later.
    """
    R, Z = np.meshgrid(r_vals, z_vals, indexing='ij')

    if geometry == "G1":
        # Concentrated near r=0.3, z=0
        T00 = np.exp(-((R - 0.3)**2 / 0.05**2 + Z**2 / 0.10**2))
    elif geometry == "G2":
        # Ring at r=0.5
        T00 = np.exp(-((R - 0.5)**2 / 0.05**2 + Z**2 / 0.10**2))
    elif geometry == "G3":
        # Two lobes along z
        T00 = (np.exp(-((R - 0.3)**2 + (Z - 0.3)**2) / 0.05**2)
             + np.exp(-((R - 0.3)**2 + (Z + 0.3)**2) / 0.05**2))
    elif geometry == "G4":
        # Broad blob
        T00 = np.exp(-(R**2 + Z**2) / 0.30**2)
    else:
        T00 = np.zeros_like(R)

    return T00 * 1e-3   # scale


def run_convergence(geometry, resolutions, r_max=1.0, z_max=1.25):
    results = []
    for Nr, Nz in resolutions:
        r_vals = np.linspace(0.0, r_max, Nr)
        z_vals = np.linspace(-z_max, z_max, Nz)

        T00 = make_test_source(r_vals, z_vals, geometry)
        h   = solve_poisson_axisym(T00, r_vals, z_vals)

        h_center = h[0, Nz // 2]
        h_max    = np.max(np.abs(h))

        results.append({
            "Nr": Nr,
            "Nz": Nz,
            "h_center": float(h_center),
            "h_max":    float(h_max),
        })
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    args = ap.parse_args()

    if args.mode == "quick":
        resolutions = [(40, 80), (60, 120), (80, 160)]
    else:
        resolutions = [(40, 80), (60, 120), (80, 160), (100, 200), (120, 240)]

    print("=" * 72)
    print("ZAMANDA YOLCULUK — ZT-006.2.21")
    print("AXISYMMETRIC POISSON SOLVER CONVERGENCE")
    print("=" * 72)
    print(f"Einstein factor 16πG/c⁴ = {EINSTEIN_FACTOR:.6e}")
    print(f"Mode: {args.mode}")
    print()

    all_results = {}
    for geom in ["G1", "G2", "G3", "G4"]:
        print(f"[{geom}]")
        res = run_convergence(geom, resolutions)
        for r in res:
            print(f"  grid=({r['Nr']:3d},{r['Nz']:3d})  "
                  f"h_center={r['h_center']:.6e}  h_max={r['h_max']:.6e}")
        for i in range(1, len(res)):
            h1 = res[i-1]["h_center"]
            h2 = res[i]["h_center"]
            if abs(h2) > 0:
                rel = abs(h2 - h1) / abs(h2)
                print(f"    rel_change({res[i-1]['Nr']} → {res[i]['Nr']}) = {rel:.4%}")
        print()
        all_results[geom] = res

    out = Path("zt006_2_21_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.2.21",
        "mode": args.mode,
        "einstein_factor": EINSTEIN_FACTOR,
        "results": all_results,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: Poisson convergence test; no GR/CTC conclusion.")


if __name__ == "__main__":
    main()