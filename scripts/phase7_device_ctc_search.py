#!/usr/bin/env python3
"""
Phase 7: Device-Generated Spacetime Search (ZT-007)

Applies the validated GR pipeline to the actual device model:
1. Compute device EM stress-energy tensor T_munu (from batch_em_source)
2. Compute retarded metric perturbation hbar_munu (from retarded_gr)
3. Search for CTCs in the resulting spacetime (from ctc_search + causal_structure)

This is the first time the full chain I -> (E,B) -> T -> g -> CTC is run
on the actual device geometry, not just benchmark metrics.

Honesty constraints:
- This is LINEARIZED gravity on Minkowski background
- No CTC or time-machine conclusion follows from this module
- The device source uses MQS approximation (known non-conservation)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from zt005.batch_em_source import compute_T_batch
from zt005.retarded_gr import retarded_hbar, C, G, KAPPA
from zt005.causal_structure import classify_curve, TIMELIKE
from zt005.benchmark_metrics import Metric
from zt005.model import ActiveCylinder, Backpack, Drive


class DeviceMetric(Metric):
    """Metric constructed from device retarded solution on Minkowski background.

    g_munu = eta_munu + h_munu
    where h_munu = hbar_munu - 1/2 eta_munu hbar (trace reverse of hbar)
    """

    name = "device_linearized"
    coords = ("t", "x", "y", "z")
    admits_ctc = False  # Linearized on Minkowski - no CTC expected

    def __init__(self, hbar_grid, obs_points, obs_times):
        """
        hbar_grid: (Nt, Nx, Ny, Nz, 4, 4) trace-reversed perturbation
        obs_points: (Nx, Ny, Nz, 3) spatial coordinates
        obs_times: (Nt,) time coordinates
        """
        self.hbar_grid = np.asarray(hbar_grid, float)
        self.obs_points = np.asarray(obs_points, float)
        self.obs_times = np.asarray(obs_times, float)
        self.Nt, self.Nx, self.Ny, self.Nz = self.hbar_grid.shape[:4]
        self.dx = float(obs_points[1, 0, 0, 0] - obs_points[0, 0, 0, 0])
        self.dy = float(obs_points[0, 1, 0, 1] - obs_points[0, 0, 0, 1])
        self.dz = float(obs_points[0, 0, 1, 2] - obs_points[0, 0, 0, 2])
        self.dt = float(obs_times[1] - obs_times[0])

    def _interp_hbar(self, x):
        """Trilinear interpolation of hbar at spacetime point x = (t, x, y, z)."""
        t, xi, yi, zi = x[0], x[1], x[2], x[3]

        # Find time index
        it = int((t - self.obs_times[0]) / self.dt)
        it = np.clip(it, 0, self.Nt - 2)

        # Find spatial indices
        ix = int((xi - self.obs_points[0, 0, 0, 0]) / self.dx)
        iy = int((yi - self.obs_points[0, 0, 0, 1]) / self.dy)
        iz = int((zi - self.obs_points[0, 0, 0, 2]) / self.dz)
        ix = np.clip(ix, 0, self.Nx - 2)
        iy = np.clip(iy, 0, self.Ny - 2)
        iz = np.clip(iz, 0, self.Nz - 2)

        # Trilinear weights
        tx = (xi - self.obs_points[ix, 0, 0, 0]) / self.dx
        ty = (yi - self.obs_points[0, iy, 0, 1]) / self.dy
        tz = (zi - self.obs_points[0, 0, iz, 2]) / self.dz

        # Interpolate in space at it and it+1
        h0 = self._trilerp(self.hbar_grid[it], ix, iy, iz, tx, ty, tz)
        h1 = self._trilerp(self.hbar_grid[it + 1], ix, iy, iz, tx, ty, tz)

        # Linear in time
        tt = (t - self.obs_times[it]) / self.dt
        return h0 + tt * (h1 - h0)

    def _trilerp(self, grid, ix, iy, iz, tx, ty, tz):
        """Trilinear interpolation on 3D grid."""
        c000 = grid[ix, iy, iz]
        c100 = grid[ix + 1, iy, iz]
        c010 = grid[ix, iy + 1, iz]
        c001 = grid[ix, iy, iz + 1]
        c110 = grid[ix + 1, iy + 1, iz]
        c101 = grid[ix + 1, iy, iz + 1]
        c011 = grid[ix, iy + 1, iz + 1]
        c111 = grid[ix + 1, iy + 1, iz + 1]

        c00 = c000 * (1 - tx) + c100 * tx
        c01 = c010 * (1 - tx) + c110 * tx
        c10 = c001 * (1 - tx) + c101 * tx
        c11 = c011 * (1 - tx) + c111 * tx

        c0 = c00 * (1 - ty) + c01 * ty
        c1 = c10 * (1 - ty) + c11 * ty

        return c0 * (1 - tz) + c1 * tz

    def g(self, x):
        """Full metric g_munu = eta_munu + h_munu at point x."""
        hbar = self._interp_hbar(x)
        # h_munu = hbar_munu - 1/2 eta_munu hbar^a_a
        tr = np.trace(hbar)
        h = hbar - 0.5 * np.diag([-1.0, 1.0, 1.0, 1.0]) * tr
        return np.diag([-1.0, 1.0, 1.0, 1.0]) + h

    def in_domain(self, x):
        return True


def build_device_source_grid(config_name="G1", n_grid=17, radius_m=0.5):
    """Build 3D grid of device stress-energy tensor at t=0."""
    bp = Backpack()
    drive = Drive(peak_current_a=100.0)
    cyl = ActiveCylinder()

    # 3D grid
    s = np.linspace(-radius_m, radius_m, n_grid)
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    # Compute T at t=0 (static snapshot)
    T = compute_T_batch(config_name, pts, bp, drive)
    return pts, T, s, X, Y, Z


def run_device_ctc_search(config_name="G1", n_grid=17, radius_m=0.5, n_time=5):
    """Run the full device CTC search pipeline."""
    print(f"=" * 80)
    print(f"PHASE 7: Device-Generated Spacetime Search")
    print(f"Config: {config_name}, Grid: {n_grid}^3, Radius: {radius_m}m")
    print(f"=" * 80)

    bp = Backpack()
    drive = Drive(peak_current_a=100.0)
    cyl = ActiveCylinder()

    # 1. Build source grid
    print("\n[1/5] Building device source grid...")
    src_pts, src_T, s, X, Y, Z = build_device_source_grid(config_name, n_grid, radius_m)
    src_weights = np.full(len(src_pts), (s[1] - s[0]) ** 3)
    print(f"    Source points: {len(src_pts)}, T shape: {src_T.shape}")

    # 2. Define observation grid (same as source for now)
    print("\n[2/5] Setting up observation grid...")
    obs_pts = src_pts.copy()
    obs_times = np.linspace(-1e-7, 1e-7, n_time)
    print(f"    Obs points: {len(obs_pts)}, Times: {len(obs_times)}")

    # 3. Compute retarded hbar
    print("\n[3/5] Computing retarded metric perturbation...")
    t0 = time.time()

    def source_fn(t_array, x_array):
        """Source function for retarded integral.
        
        t_array: (Ns,) retarded times for each source point
        x_array: (Ns, 3) source point positions (same as src_pts)
        Returns: (Ns, 4, 4) stress-energy at each source point
        """
        # For static source, T is time-independent and same for all source points
        # t_array has shape (Ns,), we need to return (Ns, 4, 4)
        # src_T already has shape (Ns, 4, 4) where Ns = len(src_pts)
        return src_T

    hbar = retarded_hbar(source_fn, obs_pts, obs_times, src_pts, src_weights)
    print(f"    hbar shape: {hbar.shape}, time: {time.time() - t0:.2f}s")
    print(f"    hbar range: [{hbar.min():.2e}, {hbar.max():.2e}]")

    # 4. Build device metric
    print("\n[4/5] Building device metric...")
    # Reshape for interpolation: (Nt, Nx, Ny, Nz, 4, 4)
    hbar_grid = hbar.reshape(n_time, n_grid, n_grid, n_grid, 4, 4)
    obs_grid = np.stack([
        np.broadcast_to(X, (n_time, n_grid, n_grid, n_grid)),
        np.broadcast_to(Y, (n_time, n_grid, n_grid, n_grid)),
        np.broadcast_to(Z, (n_time, n_grid, n_grid, n_grid)),
    ], axis=-1)

    device_metric = DeviceMetric(hbar_grid, obs_grid[0], obs_times)

    # 5. Search for CTCs (phi-circles at various radii)
    print("\n[5/5] Searching for CTCs in device spacetime...")
    from zt005.ctc_search import check_phi_circle_ctc

    # Test phi-circles at various radii in the xy-plane (z=0)
    radii = np.linspace(0.1, radius_m * 0.8, 10)
    ctc_results = []

    for r in radii:
        # Sample phi-circle at t=0, z=0
        n_phi = 65
        phi = np.linspace(0, 2 * np.pi, n_phi)
        pts = np.zeros((n_phi, 4))
        pts[:, 0] = 0.0  # t = 0
        pts[:, 1] = r * np.cos(phi)  # x
        pts[:, 2] = r * np.sin(phi)  # y
        pts[:, 3] = 0.0  # z

        # Classify using causal_structure
        rep = classify_curve(device_metric, pts, periodic_coord=1)

        # Check if closed in covering space
        d_close = pts[-1] - pts[0]
        closed = np.allclose(d_close, 0.0, atol=1e-9)

        # Check if timelike everywhere
        all_timelike = rep["overall"] == TIMELIKE

        ctc_results.append({
            "radius": float(r),
            "closed": bool(closed),
            "overall": rep["overall"],
            "is_ctc_candidate": bool(closed and all_timelike),
            "proper_time": rep["proper_time"],
        })

    # Also test coordinate circles in cylindrical-like coordinates
    # For a more proper test, we'd need cylindrical coordinates
    # But in linearized gravity on Minkowski, no CTC expected

    any_ctc = any(r["is_ctc_candidate"] for r in ctc_results)

    # Summary
    print(f"\n{'=' * 80}")
    print(f"DEVICE CTC SEARCH RESULTS")
    print(f"{'=' * 80}")
    print(f"Config: {config_name}")
    print(f"Grid: {n_grid}^3, radius: {radius_m}m")
    print(f"Max |hbar|: {np.max(np.abs(hbar)):.2e}")
    print(f"Any CTC candidate: {any_ctc}")
    print(f"\nRadius scan:")
    for r in ctc_results:
        status = "CTC!" if r["is_ctc_candidate"] else "no"
        print(f"  r={r['radius']:.3f}m: {r['overall']:12s} closed={r['closed']} -> {status}")

    # Save results
    OUT = Path("results/validated/phase7_device_ctc_search.json")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "stage": "Phase 7 - Device CTC Search",
        "config": config_name,
        "grid": {"n": n_grid, "radius_m": radius_m, "n_time": n_time},
        "hbar_stats": {
            "min": float(hbar.min()),
            "max": float(hbar.max()),
            "mean_abs": float(np.mean(np.abs(hbar))),
        },
        "ctc_search": ctc_results,
        "any_ctc_candidate": any_ctc,
        "scope_note": "Linearized gravity on Minkowski background. No CTC expected. MQS source non-conservation known.",
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nResults written to {OUT}")

    return report


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="G1", choices=["G1", "G2", "G3", "G4"])
    ap.add_argument("--n-grid", type=int, default=17)
    ap.add_argument("--radius", type=float, default=0.5)
    ap.add_argument("--n-time", type=int, default=5)
    args = ap.parse_args()

    run_device_ctc_search(args.config, args.n_grid, args.radius, args.n_time)