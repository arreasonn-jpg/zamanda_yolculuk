"""ZT-006.2.25 - retarded-solver device-source convergence.

Open item left by ZT-006.2.19/.20 (frozen kernel) and closed for the
*instantaneous* kernel by ZT-006.2.24 (a_wire_eff = 0.08 m):

    "Convergence of the retarded integral for the *device* source
    distribution"  (docs/validation.md, "Not yet validated").

This stage runs the grid-refinement study against the retarded
linearized-Einstein solver (zt005.retarded_gr.retarded_hbar) with the
actual device sources (G1-G4, external observers, same set as 2.24)
and the same effective wire smoothing a_wire_eff = 0.08 m.

Device time dependence (MQS 1 MHz drive):
    I(t) = I0 cos(w t)  =>
    A(x, t) = A1(x) cos(w t),   B(x, t) = B1(x) cos(w t),
    E(x, t) = w A1(x) sin(w t)      (E = -dA/dt)

The per-point retarded time t_ret = t - R/c produced by the solver
carries the MQS phase w t_ret ~ w*2.5 m / c ~ 5e-8 rad across the
source: formally kept, numerically far below the quadrature error.
Only the first upsampled*downsampled*up/downsampled novelty is the
retarded phase; the spatial shapes are the validated 2.24 fields.

Criterion: rel < 0.05 between consecutive grids (same as ZT-006.2.24).
"""
import argparse, json
from pathlib import Path

import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .validated_stress_energy import em_stress_energy
from .retarded_gr import retarded_hbar
from .green_convergence import focused_source_nodes
from .geometry import build_sources

OUT = Path("results/validated/zt006_2_25_results.json")
A_WIRE_EFF = 0.08      # m, same as ZT-006.2.24
REL_TOL = 0.05
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["quick", "full"], default="quick")
    a = ap.parse_args()
    cyl, bk, dv = ActiveCylinder(), Backpack(), Drive()
    omega = 2.0 * np.pi * dv.frequency_hz
    # external observers, identical to ZT-006.2.24
    obs = np.array([[0., 0., 2.], [0., 0., 3.], [0., 0., 5.],
                    [2., 0., 0.], [3., 0., 0.]])
    # observation times: a full half cycle, spread inside the cycle
    period = 1.0 / dv.frequency_hz
    times = np.array([0.0, period / 16, period / 8, period / 4,
                      3 * period / 8, period / 2])
    grids = ([(6, 12, 6), (8, 16, 8), (12, 24, 12)] if a.mode == "quick"
             else [(8, 16, 8), (12, 24, 12), (16, 32, 16), (20, 40, 20)])

    out = {"stage": "ZT-006.2.25", "mode": a.mode,
           "a_wire_effective_m": A_WIRE_EFF,
           "physical_wire_radius_m": 0.005,
           "solver": "retarded_gr.retarded_hbar",
           "obs_times_s": times.tolist(),
           "results": {}}
    all_ok = True
    for name in ["G1", "G2", "G3", "G4"]:
        srcs = build_sources(name, bk)
        fields = {}
        for g in grids:
            pts, w = focused_source_nodes(cyl.radius_m, cyl.height_m, *g)
            A1, B1 = device_field_shapes(pts, srcs, dv, A_WIRE_EFF)

            def source_fn(t_ret, _pts):
                t_ret = np.asarray(t_ret, float)[:, None]
                c = np.cos(omega * t_ret)
                s = np.sin(omega * t_ret)
                E = omega * A1 * s
                B = B1 * c
                return em_stress_energy(E, B)["T"]

            fields[str(g)] = retarded_hbar(
                source_fn, obs, times, pts, w)
        rows = []
        for i in range(1, len(grids)):
            hp = fields[str(grids[i - 1])]
            hc = fields[str(grids[i])]
            rel = float(np.linalg.norm(hc - hp)
                        / max(np.linalg.norm(hc), 1e-300))
            rows.append({"from": list(grids[i - 1]),
                         "to": list(grids[i]), "rel": rel})
            if rel >= REL_TOL:
                all_ok = False
        out["results"][name] = rows
    out["convergence_passed"] = all_ok
    out["criterion"] = f"rel < {REL_TOL} between consecutive grids"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("\nCONVERGENCE:", "PASS" if all_ok else "FAIL")


if __name__ == "__main__":
    main()
