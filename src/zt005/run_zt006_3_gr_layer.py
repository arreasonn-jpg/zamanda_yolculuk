"""ZT-006.3 — GR validation layer runner.

Runs the new physics layer end to end and writes a PASS/FAIL report:

  PHASE 3 pieces : retarded source solver, conservation diagnostic,
                   Lorenz-gauge diagnostic
  PHASE 4/6 pieces: known-metric benchmarks (Minkowski, Schwarzschild,
                   Kerr, Gödel, Tipler idealization, Morris-Thorne),
                   geodesic integrator, causal classifier, known-answer
                   CTC benchmark

No CTC or time-machine conclusion about the device is made here; the layer
validates the solver against exact known answers only.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from zt005.benchmark_metrics import all_benchmarks
from zt005.causal_structure import classify_curve
from zt005.ctc_search import (godel_ctc_benchmark, negative_benchmarks,
                              tipler_ctc_benchmark)
from zt005.gauge_checks import lorenz_residual
from zt005.geodesic_solver import integrate_geodesic
from zt005.manufactured_sources import (GaussianBump, conserved_source_fn)
from zt005.retarded_gr import C, retarded_hbar
from zt005.stress_energy_conservation import (divergence_residual,
                                              sample_source_on_grid)

OUT = Path("results/validated/zt006_3_gr_layer_results.json")


def _metric_sanity():
    probes = {
        "minkowski": np.array([0.0, 1.0, 1.0, 1.0]),
        "schwarzschild": np.array([0.0, 6.0, np.pi / 2, 0.0]),
        "kerr": np.array([0.0, 10.0, np.pi / 2, 0.0]),
        "godel": np.array([0.0, 1.2, 0.0, 0.0]),
        "tipler_cylinder": np.array([0.0, 1.5, 0.0, 0.0]),
        "morris_thorne": np.array([0.0, 0.0, np.pi / 2, 0.0]),
    }
    rows = {}
    for m in all_benchmarks():
        g = np.asarray(m.g(probes[m.name]), float)
        ev = np.linalg.eigvalsh(g)
        rows[m.name] = {
            "symmetric": bool(np.allclose(g, g.T)),
            "lorentzian_signature": bool(np.sum(ev < 0) == 1),
            "det": float(np.linalg.det(g)),
        }
    return rows


def _norm_conservation():
    from zt005.benchmark_metrics import (GodelMetric, KerrMetric,
                                         MorrisThorneMetric,
                                         SchwarzschildMetric)
    cases = {
        "schwarzschild_analytic": (SchwarzschildMetric(), [0, 10, np.pi/2, 0],
                                   [3.0, 0.0, 0.0, 0.22], True),
        "kerr_numeric": (KerrMetric(M=1.0, a=0.7), [0, 10, np.pi/3, 0],
                         [2.0, 0.1, 0.05, 0.05], False),
        "godel_numeric": (GodelMetric(), [0, 0.5, 0.3, 0.2],
                          [1.0, 0.2, 0.3, 0.1], False),
        "morris_thorne_numeric": (MorrisThorneMetric(), [0, -3, np.pi/2, 0],
                                  [1.0, 0.6, 0.0, 0.0], False),
    }
    rows = {}
    for name, (m, x0, u0_raw, analytic) in cases.items():
        x0 = np.array(x0, float)
        u0 = np.array(u0_raw, float)
        n = float(u0 @ m.g(x0) @ u0)
        u0 = u0 * np.sqrt(-1.0 / n)
        res = integrate_geodesic(m, x0, u0, tau_span=8.0, n_steps=1200,
                                 analytic=analytic)
        rows[name] = {
            "norm0": res["norm0"],
            "norm_drift_rel": res["norm_drift_rel"],
            "pass": bool(res["norm_drift_rel"] < 1e-5),
        }
    return rows


def _perihelion():
    from zt005.benchmark_metrics import SchwarzschildMetric
    M, r_p = 1.0, 40.0
    m = SchwarzschildMetric(M=M)
    x0 = np.array([0.0, r_p, np.pi / 2, 0.0])
    k = 1.08 * np.sqrt(M / r_p ** 3)
    f = 1.0 - 2.0 * M / r_p
    u_t = np.sqrt(1.0 / (f - r_p * r_p * k * k))
    u0 = np.array([u_t, 0.0, 0.0, k * u_t])
    res = integrate_geodesic(m, x0, u0, tau_span=2500.0, n_steps=50000,
                             analytic=True)
    r, phi = res["x"][:, 1], res["x"][:, 3]
    mins = [i for i in range(1, len(r) - 1)
            if r[i] < r[i - 1] and r[i] < r[i + 1]]
    r_a = float(np.max(r[:mins[0] + 1]))
    a = 0.5 * (r_p + r_a)
    e = (r_a - r_p) / (r_a + r_p)
    dphi = float(phi[mins[0]] - phi[0])
    predicted = 2 * np.pi + 6 * np.pi * M / (a * (1 - e * e))
    return {
        "semi_major_axis_M": a, "eccentricity": e,
        "measured_dphi_rad": dphi, "predicted_1pn_rad": predicted,
        "rel_err": float(abs(dphi - predicted) / predicted),
        "pass": bool(abs(dphi - predicted) / predicted < 0.05),
    }


def _conservation_gauge(mode):
    bump = GaussianBump(A=1.0, sigma_t=1e-7, sigma_x=1.0)

    n_src = 21 if mode == "full" else 19
    n_t = 13
    t = np.linspace(-3e-7, 3e-7, 25)
    x = np.linspace(-3.0, 3.0, 25)
    T = sample_source_on_grid(conserved_source_fn(bump), t, x, x, x)
    cons = divergence_residual(T, t[1] - t[0], x[1] - x[0])

    s = np.linspace(-3.0, 3.0, n_src)
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    w = np.full(len(pts), (s[1] - s[0]) ** 3)
    o = np.linspace(-0.6, 0.6, 7)
    OX, OY, OZ = np.meshgrid(o, o, o, indexing="ij")
    obs = np.stack([OX.ravel(), OY.ravel(), OZ.ravel()], axis=1)
    times = np.linspace(-6e-8, 6e-8, n_t)
    h = retarded_hbar(conserved_source_fn(bump), obs, times, pts, w)
    grid = h.reshape(len(times), 7, 7, 7, 4, 4)
    lz = lorenz_residual(grid, times[1] - times[0], o[1] - o[0])
    L = min(C * bump.st, bump.sx)
    return {
        "conservation_max_rel": cons["max_rel"],
        "conservation_pass": bool(cons["max_rel"] < 2e-2),
        "lorenz_scaled_residual": float(lz["max_abs"] * L / lz["scale"]),
        "lorenz_pass": bool(lz["max_abs"] * L / lz["scale"] < 0.15),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["quick", "full"], default="quick")
    a = ap.parse_args()

    t0 = time.time()
    report = {
        "stage": "ZT-006.3",
        "layer": "GR validation layer (phases 3/4/6 groundwork)",
        "mode": a.mode,
        "scope_note": ("Known-answer validation only. No CTC or device "
                       "conclusion is produced."),
    }

    print("=" * 76)
    print("ZT-006.3 — GR VALIDATION LAYER")
    print("=" * 76)

    report["metric_sanity"] = _metric_sanity()
    print("[1/5] metric sanity:",
          {k: v["lorentzian_signature"] for k, v
           in report["metric_sanity"].items()})

    report["geodesic_norm_conservation"] = _norm_conservation()
    print("[2/5] norm drift:",
          {k: v["norm_drift_rel"] for k, v
           in report["geodesic_norm_conservation"].items()})

    report["schwarzschild_perihelion"] = _perihelion()
    print("[3/5] perihelion rel err:",
          report["schwarzschild_perihelion"]["rel_err"])

    report["ctc_known_answer"] = {
        "godel": godel_ctc_benchmark(n_points=97),
        "tipler": tipler_ctc_benchmark(n_points=97),
        "negative": negative_benchmarks(n_points=49),
    }
    for key in ("godel", "tipler"):
        d = report["ctc_known_answer"][key]
        for part in ("inside", "outside"):
            keep = ("overall", "proper_time", "closed_covering_space",
                    "is_ctc_candidate")
            d[part] = {k: v for k, v in d[part].items() if k in keep}
    print("[4/5] CTC known-answer:",
          report["ctc_known_answer"]["godel"]["agrees_with_analytic"],
          report["ctc_known_answer"]["tipler"]["agrees_with_analytic"],
          report["ctc_known_answer"]["negative"]["all_clean"])

    report["conservation_and_gauge"] = _conservation_gauge(a.mode)
    print("[5/5] conservation/gauge:", report["conservation_and_gauge"])

    checks = (
        [v["symmetric"] and v["lorentzian_signature"]
         for v in report["metric_sanity"].values()]
        + [v["pass"] for v in report["geodesic_norm_conservation"].values()]
        + [report["schwarzschild_perihelion"]["pass"]]
        + [report["ctc_known_answer"]["godel"]["agrees_with_analytic"],
           report["ctc_known_answer"]["tipler"]["agrees_with_analytic"],
           report["ctc_known_answer"]["negative"]["all_clean"]]
        + [report["conservation_and_gauge"]["conservation_pass"],
           report["conservation_and_gauge"]["lorenz_pass"]]
    )
    report["overall"] = {
        "n_checks": len(checks),
        "n_pass": int(sum(checks)),
        "status": "PASS" if all(checks) else "FAIL",
        "runtime_s": round(time.time() - t0, 2),
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, default=float),
                   encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("OVERALL:", report["overall"])
    print("No CTC conclusion about the device is made at this stage.")


if __name__ == "__main__":
    main()
