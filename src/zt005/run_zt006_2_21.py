"""ZT-006.2.21 — Current Normalization, Dual Independent Green Solvers, and Convergence Gate.

Milestone runner closing Section I requirements:
1. Current normalization audit: I_parent, I_total, I_peak, I_rms across quadrature resolutions.
2. Dual independent Green solvers: Solver A (Gauss-cylindrical) vs Solver B (Cartesian voxel with exact boundary flux).
3. Convergence gate: grid series refinement with final delta < 1% on h00, |h0i|, and ||hij||.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .current_normalization import multi_resolution_current_audit
from .green_convergence import run_convergence_study
from .green_solvers import (
    solver_a_tensor_product,
    solver_b_adaptive_voxel,
    compare_green_solvers,
)
from .gauss_cyl_quadrature import gauss_cylindrical_nodes
from .manufactured_sources import GaussianBump, conserved_T

OUT = Path("results/validated/zt006_2_21_results.json")


def main():
    parser = argparse.ArgumentParser(description="ZT-006.2.21 Runner")
    parser.add_argument("--mode", choices=["quick", "full"], default="quick")
    args = parser.parse_args()

    t_start = time.time()
    print("=" * 88)
    print("ZAMANDA YOLCULUK — ZT-006.2.21")
    print("CURRENT NORMALIZATION, DUAL GREEN SOLVERS, & CONVERGENCE GATE")
    print(f"Mode: {args.mode}")
    print("=" * 88)

    cyl = ActiveCylinder()
    bp = Backpack()
    drive = Drive(peak_current_a=100.0)

    # 1. Current Normalization Audit
    print("\n--- [GATE 1/3] AKIM NORMALIZASYONU AUDIT ---")
    current_audits = {}
    current_all_pass = True
    geometries = ["G1", "G2", "G3", "G4"]
    rules = [(1, 6), (2, 8), (3, 10)] if args.mode == "quick" else [(1, 6), (2, 8), (3, 10), (4, 12), (5, 16)]

    for geom in geometries:
        audit = multi_resolution_current_audit(geom, bp, radius_m=0.0025, rules=rules, i_peak=drive.peak_current_a)
        current_audits[geom] = audit
        if not audit["all_pass"]:
            current_all_pass = False
        print(f"  {geom}: I_parent={audit['current_spec']['I_parent']} A, "
              f"I_total={audit['current_spec']['I_total']} A, "
              f"max disc={audit['max_resolution_discrepancy']:.2e} -> "
              f"{'PASS' if audit['all_pass'] else 'FAIL'}")

    # 2. Dual Independent Green Solver Validation
    print("\n--- [GATE 2/3] BAĞIMSIZ İKİNCİ GREEN SOLVER DOĞRULAMASI ---")
    val_obs = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    bump = GaussianBump(A=1.0, sigma_t=1e-3, sigma_x=0.5)

    n_a = 28 if args.mode == "quick" else 36
    pts_a, w_a = gauss_cylindrical_nodes(cyl.radius_m, cyl.height_m, n_a, 2 * n_a, n_a)
    T_a = conserved_T(bump, np.zeros(len(pts_a)), pts_a)
    h_solver_a = solver_a_tensor_product(val_obs, pts_a, w_a, T_a)

    def eval_src(pts):
        return conserved_T(bump, np.zeros(len(pts)), pts)

    n_b = 36 if args.mode == "quick" else 48
    h_solver_b = solver_b_adaptive_voxel(
        val_obs, cyl.radius_m, cyl.height_m, eval_src, nx=n_b, ny=n_b, nz=n_b
    )
    eps_tol = 0.01
    solver_comparison = compare_green_solvers(h_solver_a, h_solver_b, eps_tol=eps_tol)
    print(f"  Solver A (Tensor Gauss) vs Solver B (Adaptive Voxel):")
    print(f"    Tensor Rel Error: {solver_comparison['tensor_rel_err']:.4e} (threshold: {solver_comparison['eps_tol']})")
    print(f"    h00 Rel Error   : {solver_comparison['h00_rel_err']:.4e}")
    print(f"    hij Rel Error   : {solver_comparison['hij_rel_err']:.4e}")
    print(f"    h0i Abs Error   : {solver_comparison['h0i_abs_err']:.4e}")
    print(f"    Status          : {'PASS' if solver_comparison['pass'] else 'FAIL'}")

    # 3. Green Convergence Gate (Grid Series N=12, 16, 20, 28, 36)
    print("\n--- [GATE 3/3] GREEN-INTEGRAL YAKINSAMASI (GRID SERİSİ) ---")
    grid_ns = [12, 16, 20] if args.mode == "quick" else [12, 16, 20, 28, 36]
    conv_study = run_convergence_study(
        name="G1",
        backpack=bp,
        drive=drive,
        cyl=cyl,
        grid_ns=grid_ns,
        tol_h00=0.01,
        tol_hij=0.01,
        source_type="benchmark",
        validate_with_solver_b=False,
    )
    acc = conv_study["acceptance"]
    print(f"  Grid series evaluated: {grid_ns}")
    for step in conv_study["consecutive_steps"]:
        print(f"    Step N={step['from_N']} -> N={step['to_N']}: "
              f"max_h00_rel={step['max_h00_rel']:.4e}, "
              f"max_hij_rel={step['max_hij_rel']:.4e}, "
              f"max_h0i_abs={step['max_h0i_abs']:.4e}")

    print(f"  Final step acceptance:")
    print(f"    h00_rel: {acc['h00_rel']:.4e} (threshold: {acc['h00_threshold']}) -> {'PASS' if acc['pass_h00'] else 'FAIL'}")
    print(f"    hij_rel: {acc['hij_rel']:.4e} (threshold: {acc['hij_threshold']}) -> {'PASS' if acc['pass_hij'] else 'FAIL'}")
    print(f"    h0i_abs: {acc['h0i_abs']:.4e} -> {'PASS' if acc['pass_h0i'] else 'FAIL'}")

    overall_pass = bool(current_all_pass and solver_comparison["pass"] and acc["overall_pass"])
    runtime = time.time() - t_start

    print("\n" + "=" * 88)
    print(f"ZT-006.2.21 SUMMARY: {'ALL GATES PASSED' if overall_pass else 'SOME GATES FAILED'}")
    print(f"Runtime: {runtime:.2f} s")
    print("=" * 88)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "stage": "ZT-006.2.21",
        "mode": args.mode,
        "runtime_s": runtime,
        "overall_pass": overall_pass,
        "gate_1_current_normalization": {
            "pass": current_all_pass,
            "audits": current_audits,
        },
        "gate_2_dual_green_solvers": {
            "pass": solver_comparison["pass"],
            "comparison": solver_comparison,
        },
        "gate_3_green_convergence": {
            "pass": acc["overall_pass"],
            "convergence_study": conv_study,
        },
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Report successfully saved to {OUT}")


if __name__ == "__main__":
    main()
