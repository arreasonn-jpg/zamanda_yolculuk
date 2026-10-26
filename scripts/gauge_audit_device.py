#!/usr/bin/env python3
"""
Gauge audit of the device pipeline: d_mu hbar^munu for the device retarded solution.

This is the critical open item from the roadmap (Phase 3 closure):
"The exact retarded solution of a CONSERVED source satisfies
partial_mu hbar^{mu nu} = 0 identically."

We compute the Lorenz-gauge residual for the actual device source
using the retarded solver (zt005.retarded_gr.retarded_hbar) and
the gauge diagnostic (zt005.gauge_checks.lorenz_residual).
"""

import argparse
import json
from pathlib import Path

import numpy as np

from zt005.model import ActiveCylinder, Backpack, Drive
from zt005.geometry import build_sources
from zt005.green_convergence import focused_source_nodes
from zt005.retarded_gr import retarded_hbar
from zt005.validated_stress_energy import em_stress_energy
from zt005.gauge_checks import lorenz_residual


# Same effective wire smoothing as ZT-006.2.24/2.25
A_WIRE_EFF = 0.08
CHUNK = 8000


def device_field_shapes(points, sources, drive: Drive, a_wire: float):
    """A1, B1 at peak current (instantaneous MQS shapes), no time factor."""
    p = np.asarray(points, float)
    N = len(p)
    A = np.zeros_like(p)
    B = np.zeros_like(p)
    for pts_s, seg in sources:
        mid = pts_s + 0.5 * seg
        dlI = seg * drive.peak_current_a
        for i in range(0, N, CHUNK):
            ch = p[i:i + CHUNK]
            r = ch[:, None, :] - mid[None, :, :]
            R_reg = np.sqrt(np.sum(r * r, axis=2) + a_wire ** 2)
            A[i:i + CHUNK] += (1e-7) * np.sum(
                dlI[None, :, :] / R_reg[:, :, None], axis=1)
            B[i:i + CHUNK] += (1e-7) * np.sum(
                np.cross(dlI[None, :, :], r, axis=2)
                / R_reg[:, :, None] ** 3, axis=1)
    return A, B


def build_spacetime_grid(obs_center, obs_radius, n_obs, n_time, period):
    """Build a 4D spacetime grid for Lorenz residual computation."""
    # Spatial grid: cubic around obs_center
    o = np.linspace(-obs_radius, obs_radius, n_obs)
    OX, OY, OZ = np.meshgrid(o, o, o, indexing="ij")
    obs_pts = np.stack([OX.ravel(), OY.ravel(), OZ.ravel()], axis=1)
    obs_pts += np.asarray(obs_center, float)

    # Time grid: centered around t=0, spanning a few periods
    times = np.linspace(-1.5 * period, 1.5 * period, n_time)

    return obs_pts, times, o[1] - o[0], times[1] - times[0]


def run_gauge_audit(config_name, grid_size, mode="quick"):
    """Run gauge audit for a specific device configuration."""
    cyl, bk, dv = ActiveCylinder(), Backpack(), Drive()
    omega = 2.0 * np.pi * dv.frequency_hz
    period = 1.0 / dv.frequency_hz

    # Source grid
    if mode == "quick":
        src_grids = [(6, 12, 6), (8, 16, 8)]
    else:
        src_grids = [(8, 16, 8), (12, 24, 12), (16, 32, 16)]

    src_grid = src_grids[grid_size]
    pts, w = focused_source_nodes(cyl.radius_m, cyl.height_m, *src_grid)
    srcs = build_sources(config_name, bk)
    A1, B1 = device_field_shapes(pts, srcs, dv, A_WIRE_EFF)

    def source_fn(t_ret, _pts):
        t_ret = np.asarray(t_ret, float)[:, None]
        c = np.cos(omega * t_ret)
        s = np.sin(omega * t_ret)
        E = omega * A1 * s
        B = B1 * c
        return em_stress_energy(E, B)["T"]

    # Observation grid for gauge check
    # Use a region around the device (z-axis, near the cylinder)
    obs_center = [0.0, 0.0, 0.0]
    obs_radius = cyl.radius_m * 1.5  # 1.5x cylinder radius
    n_obs = 9  # 9^3 = 729 points
    n_time = 13

    obs_pts, times, dx, dt = build_spacetime_grid(
        obs_center, obs_radius, n_obs, n_time, period)

    print(f"  Computing retarded hbar on {len(times)}x{len(obs_pts)} grid...")
    h = retarded_hbar(source_fn, obs_pts, times, pts, w)

    # Reshape to (Nt, Nx, Ny, Nz, 4, 4)
    nt, no = len(times), n_obs
    h_grid = h.reshape(nt, no, no, no, 4, 4)

    print(f"  Computing Lorenz residual...")
    rep = lorenz_residual(h_grid, dt, dx)

    L = min(299792458.0 * 1e-6, cyl.radius_m)  # c*sigma_t or sigma_x
    scaled = rep["max_abs"] * L / rep["scale"]

    return {
        "config": config_name,
        "src_grid": src_grid,
        "obs_grid": {"center": obs_center, "radius": obs_radius,
                     "n_obs": n_obs, "n_time": n_time},
        "dx_m": dx,
        "dt_s": dt,
        "max_abs_residual": rep["max_abs"],
        "max_rel_residual": rep["max_rel"],
        "per_nu_abs": rep["per_nu_abs"],
        "hbar_scale": rep["scale"],
        "natural_length_scale": L,
        "scaled_residual": scaled,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Gauge audit: d_mu hbar^munu for device retarded solution")
    ap.add_argument("--mode", choices=["quick", "full"], default="quick")
    ap.add_argument("--configs", nargs="+",
                    default=["G1", "G2", "G3", "G4"],
                    help="Device configurations to test")
    ap.add_argument("--grid-size", type=int, default=0,
                    help="Source grid index (0=coarsest)")
    ap.add_argument("--output", type=Path,
                    default=Path("results/validated/gauge_audit_device.json"))
    a = ap.parse_args()

    print(f"Gauge audit of device pipeline (mode={a.mode})")
    print(f"Configs: {a.configs}")
    print(f"Source grid index: {a.grid_size}")
    print()

    results = {"stage": "gauge_audit_device", "mode": a.mode,
               "a_wire_eff_m": A_WIRE_EFF,
               "results": {}}

    for name in a.configs:
        print(f"[{name}] Running gauge audit...")
        try:
            res = run_gauge_audit(name, a.grid_size, a.mode)
            results["results"][name] = res
            print(f"  max|d_mu hbar^munu| = {res['max_abs_residual']:.3e}")
            print(f"  max|d_mu hbar^munu|/|hbar| = {res['max_rel_residual']:.3e}")
            print(f"  scaled_residual = {res['scaled_residual']:.3e}")
            print(f"  per-nu: {[f'{x:.3e}' for x in res['per_nu_abs']]}")
        except Exception as e:
            print(f"  ERROR: {e}")
            results["results"][name] = {"error": str(e)}
        print()

    # Summary
    print("=" * 60)
    print("SUMMARY")
    for name, res in results["results"].items():
        if "error" in res:
            print(f"  {name}: ERROR - {res['error']}")
        else:
            status = "PASS" if res["scaled_residual"] < 0.15 else "FAIL"
            print(f"  {name}: scaled_residual={res['scaled_residual']:.3e} [{status}]")

    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(results, indent=2))
    print(f"\nResults written to {a.output}")


if __name__ == "__main__":
    main()