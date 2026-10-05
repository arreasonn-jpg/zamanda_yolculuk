"""Auto-patch script for ZT-006.2.22. Run once: python _patch_zt006_2_22.py"""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_FILE  = PROJECT_ROOT / "src" / "zt005" / "run_zt006_2_22.py"
TEST_FILE = PROJECT_ROOT / "tests" / "test_zt006_2_22.py"

SRC_CODE = r'''
"""
ZT-006.2.22 - Full chain: 3D finite-wire source -> axisym T00 -> Poisson -> h_munu.
"""
import argparse
import json
from pathlib import Path
import numpy as np

from zt005.model import Backpack
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.poisson_solver import solve_poisson_axisym, EINSTEIN_FACTOR


def _flatten_sources(sources_grouped):
    """
    Extract all 3D wire points from build_finite_sources_grouped output.

    Actual structure:
      list of N groups (one per wire)
      each group: np.ndarray shape (n_fil, 2, n_seg, 3)
      axis 1 = [base_path, offset_or_direction]
      we use axis1=0 (base_path) and ignore offset info
    """
    pts_list = []
    for group in sources_grouped:
        arr = np.asarray(group, dtype=float)
        if arr.size == 0:
            continue
        if arr.ndim == 4 and arr.shape[-1] == 3:
            base = arr[:, 0, :, :]
            pts_list.append(base.reshape(-1, 3))
        elif arr.ndim == 3 and arr.shape[-1] == 3:
            pts_list.append(arr.reshape(-1, 3))
        elif arr.ndim == 2 and arr.shape[-1] == 3:
            pts_list.append(arr)

    if not pts_list:
        return np.zeros((0, 3)), np.zeros(0)
    all_pts = np.vstack(pts_list)
    currents = np.ones(len(all_pts))
    return all_pts, currents


def build_T00_axisym_smoothed(geometry, r_vals, z_vals, radius_m=0.005,
                              sigma_smooth=0.02, backpack=None,
                              current_scale=1.0):
    """
    Build 2D axisymmetric T00(r,z) from 3D finite-wire source via
    Gaussian smoothing. Backpack center is treated as axisym origin.
    """
    if backpack is None:
        backpack = Backpack()

    sources = build_finite_sources_grouped(geometry, backpack, radius_m=radius_m)
    pts, currents = _flatten_sources(sources)

    Nr, Nz = len(r_vals), len(z_vals)
    T00 = np.zeros((Nr, Nz))
    if len(pts) == 0:
        return T00

    y_center = float(getattr(backpack, "center_y_m", 0.0))
    r_pts = np.sqrt(pts[:, 0]**2 + (pts[:, 1] - y_center)**2)
    z_pts = pts[:, 2]
    weights = (currents * current_scale)**2

    two_sigma2 = 2.0 * sigma_smooth**2
    norm = 2.0 * np.pi * sigma_smooth**2
    R, Z = np.meshgrid(r_vals, z_vals, indexing="ij")

    chunk = 2000
    for s in range(0, len(pts), chunk):
        e = min(s + chunk, len(pts))
        rp = r_pts[s:e]; zp = z_pts[s:e]; wp = weights[s:e]
        d2 = (R[None] - rp[:, None, None])**2 + (Z[None] - zp[:, None, None])**2
        k = np.exp(-d2 / two_sigma2)
        T00 += np.einsum("i,ijk->jk", wp, k) / norm

    return T00


def run_convergence(geometry, resolutions, r_max=1.0, z_max=1.25, sigma=0.02):
    results = []
    for Nr, Nz in resolutions:
        r_vals = np.linspace(0.0, r_max, Nr)
        z_vals = np.linspace(-z_max, z_max, Nz)
        T00 = build_T00_axisym_smoothed(geometry, r_vals, z_vals, sigma_smooth=sigma)
        h = solve_poisson_axisym(T00, r_vals, z_vals)
        results.append({
            "Nr": Nr, "Nz": Nz,
            "T00_max": float(np.max(T00)),
            "h_center": float(h[0, Nz // 2]),
            "h_max": float(np.max(np.abs(h))),
        })
    return results


def _print_diagnostics(geometry):
    b = Backpack()
    sources = build_finite_sources_grouped(geometry, b, radius_m=0.005)
    pts, currents = _flatten_sources(sources)
    yc = float(getattr(b, "center_y_m", 0.0))
    print(f"  [{geometry}] N_groups={len(sources)}  N_points={len(pts)}")
    if len(pts) > 0:
        r = np.sqrt(pts[:, 0]**2 + (pts[:, 1] - yc)**2)
        print(f"    x: [{pts[:,0].min():.4f}, {pts[:,0].max():.4f}]")
        print(f"    y: [{pts[:,1].min():.4f}, {pts[:,1].max():.4f}]")
        print(f"    z: [{pts[:,2].min():.4f}, {pts[:,2].max():.4f}]")
        print(f"    y_center={yc:.4f}  shifted r: [{r.min():.4f}, {r.max():.4f}]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--sigma", type=float, default=0.02)
    args = ap.parse_args()

    if args.mode == "quick":
        resolutions = [(40, 80), (60, 120), (80, 160)]
    else:
        resolutions = [(40, 80), (60, 120), (80, 160), (100, 200)]

    print("=" * 72)
    print("ZAMANDA YOLCULUK - ZT-006.2.22")
    print("FULL CHAIN: 3D FINITE-WIRE -> AXISYM T00 -> POISSON -> h_munu")
    print("=" * 72)
    print(f"Einstein factor 16piG/c^4 = {EINSTEIN_FACTOR:.6e}")
    print(f"Gaussian sigma = {args.sigma} m")
    print(f"Mode: {args.mode}")
    print()

    print("Source diagnostics:")
    for geom in ["G1", "G2", "G3", "G4"]:
        try:
            _print_diagnostics(geom)
        except Exception as e:
            print(f"  [{geom}] FAILED: {e}")
    print()

    all_results = {}
    for geom in ["G1", "G2", "G3", "G4"]:
        print(f"[{geom}]")
        try:
            res = run_convergence(geom, resolutions, sigma=args.sigma)
            for r in res:
                print(f"  grid=({r['Nr']:3d},{r['Nz']:3d})  "
                      f"T00_max={r['T00_max']:.4e}  "
                      f"h_center={r['h_center']:.6e}  "
                      f"h_max={r['h_max']:.6e}")
            for i in range(1, len(res)):
                h1 = res[i-1]["h_center"]
                h2 = res[i]["h_center"]
                if abs(h2) > 0:
                    rel = abs(h2 - h1) / abs(h2)
                    print(f"    rel_change({res[i-1]['Nr']} -> {res[i]['Nr']}) = {rel:.4%}")
            all_results[geom] = res
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}")
            all_results[geom] = []
        print()

    out = Path("zt006_2_22_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.2.22",
        "mode": args.mode,
        "sigma": args.sigma,
        "einstein_factor": EINSTEIN_FACTOR,
        "results": all_results,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: full-chain Poisson convergence; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TEST_CODE = r'''
"""ZT-006.2.22 unit tests."""
import numpy as np
from zt005.run_zt006_2_22 import (
    _flatten_sources,
    build_T00_axisym_smoothed,
)
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.poisson_solver import solve_poisson_axisym
from zt005.model import Backpack


def test_flatten_sources_returns_nonempty_points():
    b = Backpack()
    sources = build_finite_sources_grouped("G1", b, radius_m=0.005)
    pts, currents = _flatten_sources(sources)
    assert pts.ndim == 2
    assert pts.shape[1] == 3
    assert len(pts) > 100, f"Expected >100 points, got {len(pts)}"
    assert len(currents) == len(pts)
    assert np.all(currents > 0)


def test_T00_nonzero_for_all_geometries():
    b = Backpack()
    r = np.linspace(0, 1, 30)
    z = np.linspace(-1.25, 1.25, 60)
    for geom in ["G1", "G2", "G3", "G4"]:
        T00 = build_T00_axisym_smoothed(geom, r, z, backpack=b)
        assert np.max(T00) > 0, f"{geom} produced zero T00"


def test_T00_shapes_match_grid():
    r = np.linspace(0, 1, 25)
    z = np.linspace(-1.25, 1.25, 50)
    T00 = build_T00_axisym_smoothed("G1", r, z)
    assert T00.shape == (25, 50)


def test_poisson_from_real_source_is_finite():
    r = np.linspace(0, 1, 30)
    z = np.linspace(-1.25, 1.25, 60)
    T00 = build_T00_axisym_smoothed("G1", r, z)
    h = solve_poisson_axisym(T00, r, z)
    assert np.all(np.isfinite(h))
'''

if __name__ == "__main__":
    SRC_FILE.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC_FILE}")
    TEST_FILE.write_text(TEST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TEST_FILE}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt006_2_22.py -v")
    print("  python -m zt005.run_zt006_2_22 --mode quick")