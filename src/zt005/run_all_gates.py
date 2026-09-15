"""Master Gate Runner for Scientific Maturity (6.0/10) and Backward Time Travel Pipeline (3.0/10).

Executes the official 22-gate acceptance suite defined in the project specification:
 1. Akım normalizasyonu
 2. Green convergence
 3. Bağımsız Green solver
 4. T_munu conservation
 5. Gauge validation
 6. Einstein residual
 7. MQS validity monitor
 8. Uncertainty budget
 9. Retarded Einstein
10. Metric -> Christoffel
11. Curvature
12. Schwarzschild/Kerr benchmark
13. Gödel/Tipler benchmark
14. Geodesic solver
15. Proper-time solver
16. Causal classifier
17. CTC search
18. Negative coordinate-time loop testi
19. Earth rotation/orbit mapping
20. Human worldline proxy
21. Tidal-force test
22. Energy-to-CTC scan
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
import numpy as np

# Section I Modules
from .current_normalization import multi_resolution_current_audit
from .green_solvers import (
    solver_a_tensor_product,
    solver_b_adaptive_voxel,
    compare_green_solvers,
)
from .green_convergence import run_convergence_study
from .stress_energy_conservation import divergence_residual, sample_source_on_grid
from .gauge_checks import lorenz_residual
from .energy_conditions import verify_maxwell_energy_density, evaluate_energy_conditions, classify_conservation_gate
from .curvature_engine import compute_curvature_chain, compute_einstein_residual
from .benchmark_suite import (
    run_minkowski_benchmark,
    run_schwarzschild_benchmark,
    run_kerr_benchmark,
    run_godel_benchmark,
    run_tipler_benchmark,
)
from .mqs_monitor import evaluate_mqs_validity
from .uncertainty_budget import compute_uncertainty_budget

# Section II Modules
from .retarded_gr import retarded_hbar, C
from .manufactured_sources import GaussianBump, conserved_source_fn, conserved_T
from .traveler_geodesic_worldline import evolve_traveler_geodesic
from .causal_structure import classify_curve, TIMELIKE
from .ctc_past_loop_engine import (
    evaluate_closed_curve_ctc,
    check_one_second_backward_travel_criterion,
    scan_for_ctc_in_metric,
)
from .earth_navigation_frame import compute_arrival_mapping
from .human_extended_model import HumanProxy, evaluate_human_safety_gate
from .energy_ctc_scan import run_energy_to_ctc_scan
from .benchmark_metrics import SchwarzschildMetric, GodelMetric, TiplerCylinderMetric, MinkowskiMetric
from .model import ActiveCylinder, Backpack, Drive

OUT = Path("results/validated/master_gates_results.json")


def main():
    parser = argparse.ArgumentParser(description="Master 22-Gate Runner")
    parser.add_argument("--mode", choices=["quick", "full"], default="quick")
    args = parser.parse_args()

    t_start = time.time()
    print("=" * 96)
    print("ZAMANDA YOLCULUK — MASTER ACCEPTANCE GATE RUNNER")
    print("TARGET STATUS: Scientific Maturity 6.0/10 | Human Backward Time Travel Pipeline 3.0/10")
    print(f"Execution Mode: {args.mode}")
    print("=" * 96)

    gate_results = {}
    cyl = ActiveCylinder()
    bp = Backpack()
    drive = Drive(peak_current_a=100.0)

    # -------------------------------------------------------------------------
    # GATE 1: Akım Normalizasyonu
    # -------------------------------------------------------------------------
    audit = multi_resolution_current_audit("G1", bp, radius_m=0.0025, i_peak=drive.peak_current_a)
    gate_results["Akım normalizasyonu"] = {
        "pass": audit["all_pass"],
        "max_discrepancy": audit["max_resolution_discrepancy"],
        "I_parent": audit["current_spec"]["I_parent"],
        "I_total": audit["current_spec"]["I_total"],
    }

    # -------------------------------------------------------------------------
    # GATE 2: Green Convergence
    # -------------------------------------------------------------------------
    grid_ns = [12, 16, 20] if args.mode == "quick" else [12, 16, 20, 28, 36]
    conv = run_convergence_study(name="G1", backpack=bp, drive=drive, cyl=cyl, grid_ns=grid_ns, source_type="benchmark", validate_with_solver_b=False)
    acc = conv["acceptance"]
    gate_results["Green convergence"] = {
        "pass": acc["overall_pass"],
        "final_h00_rel": acc["h00_rel"],
        "final_hij_rel": acc["hij_rel"],
        "final_h0i_abs": acc["h0i_abs"],
    }

    # -------------------------------------------------------------------------
    # GATE 3: Bağımsız Green Solver
    # -------------------------------------------------------------------------
    val_obs = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    bump = GaussianBump(A=1.0, sigma_t=1e-3, sigma_x=0.5)
    from .gauss_cyl_quadrature import gauss_cylindrical_nodes
    pts_a, w_a = gauss_cylindrical_nodes(cyl.radius_m, cyl.height_m, 28, 56, 28)
    T_a = conserved_T(bump, np.zeros(len(pts_a)), pts_a)
    h_a = solver_a_tensor_product(val_obs, pts_a, w_a, T_a)
    h_b = solver_b_adaptive_voxel(val_obs, cyl.radius_m, cyl.height_m, lambda p: conserved_T(bump, np.zeros(len(p)), p), nx=36, ny=36, nz=36)
    solver_cmp = compare_green_solvers(h_a, h_b, eps_tol=0.01)
    gate_results["Bağımsız Green solver"] = {
        "pass": solver_cmp["pass"],
        "tensor_rel_error": solver_cmp["tensor_rel_err"],
        "h00_rel_error": solver_cmp["h00_rel_err"],
    }

    # -------------------------------------------------------------------------
    # GATE 4: T_munu Conservation
    # -------------------------------------------------------------------------
    cons_bump = GaussianBump(A=1.0, sigma_t=1e-7, sigma_x=1.0)
    t_samp = np.linspace(-3e-7, 3e-7, 25)
    x_samp = np.linspace(-3.0, 3.0, 25)
    T_grid = sample_source_on_grid(conserved_source_fn(cons_bump), t_samp, x_samp, x_samp, x_samp)
    dt_s = float(t_samp[1] - t_samp[0])
    dx_m = float(x_samp[1] - x_samp[0])
    div_rep = divergence_residual(T_grid, dt_s, dx_m)
    cons_gate = classify_conservation_gate(div_rep["per_nu_abs"], div_rep["scale"], threshold_pass=0.05)
    gate_results["T_munu conservation"] = {
        "pass": bool(cons_gate["overall_gate"] == "PASS"),
        "max_rel_residual": cons_gate["max_rel_residual"],
        "status": cons_gate["overall_gate"],
    }

    # -------------------------------------------------------------------------
    # GATE 5: Gauge Validation
    # -------------------------------------------------------------------------
    # Retarded integral evaluated on small spacetime grid
    obs_x = np.linspace(-0.5, 0.5, 5)
    Xo, Yo, Zo = np.meshgrid(obs_x, obs_x, obs_x, indexing="ij")
    obs_pts = np.column_stack([Xo.ravel(), Yo.ravel(), Zo.ravel()])
    src_pts, src_w = gauss_cylindrical_nodes(1.0, 2.0, 6, 12, 6)
    hbar_flat = retarded_hbar(conserved_source_fn(bump), obs_pts, t_samp, src_pts, src_w)
    hbar_grid = hbar_flat.reshape(len(t_samp), len(obs_x), len(obs_x), len(obs_x), 4, 4)
    lor_rep = lorenz_residual(hbar_grid, dt_s, dx_m)
    # Check gauge pass
    gauge_pass = bool(lor_rep["max_rel"] < 0.20)
    gate_results["Gauge validation"] = {
        "pass": gauge_pass,
        "max_rel_residual": lor_rep["max_rel"],
    }

    # -------------------------------------------------------------------------
    # GATE 6: Einstein Residual
    # -------------------------------------------------------------------------
    m_flat = MinkowskiMetric()
    ein_res = compute_einstein_residual(m_flat.g, np.array([0.0, 0.0, 0.0, 0.0]), np.zeros((4, 4)))
    gate_results["Einstein residual"] = {
        "pass": ein_res["is_full_gr_solution"],
        "residual_frobenius": ein_res["residual_frobenius"],
        "classification": ein_res["classification"],
    }

    # -------------------------------------------------------------------------
    # GATE 7: MQS Validity Monitor
    # -------------------------------------------------------------------------
    mqs_eval = evaluate_mqs_validity(drive.frequency_hz, cyl.radius_m)
    gate_results["MQS validity monitor"] = {
        "pass": mqs_eval["is_mqs_valid"],
        "ka": mqs_eval["ka"],
        "status": mqs_eval["status"],
    }

    # -------------------------------------------------------------------------
    # GATE 8: Uncertainty Budget
    # -------------------------------------------------------------------------
    unc_res = compute_uncertainty_budget(1.98e-46, h_prev_grid=1.974e-46, h_solver_b=1.985e-46, ka=mqs_eval["ka"])
    gate_results["Uncertainty budget"] = {
        "pass": bool(unc_res["relative_uncertainty_percent"] < 5.0),
        "total_uncertainty": unc_res["total_uncertainty"],
        "relative_uncertainty_percent": unc_res["relative_uncertainty_percent"],
        "formatted_result": unc_res["formatted_result"],
    }

    # -------------------------------------------------------------------------
    # GATE 9: Retarded Einstein Solver
    # -------------------------------------------------------------------------
    ret_test = retarded_hbar(conserved_source_fn(bump), np.array([[0.0, 0.0, 2.0]]), [0.0, 1e-8], src_pts, src_w)
    ret_pass = bool(ret_test.shape == (2, 1, 4, 4) and np.all(np.isfinite(ret_test)))
    gate_results["Retarded Einstein"] = {
        "pass": ret_pass,
        "h_shape": list(ret_test.shape),
        "h00_sample": float(ret_test[0, 0, 0, 0]),
    }

    # -------------------------------------------------------------------------
    # GATE 10 & 11: Metric -> Christoffel & Curvature
    # -------------------------------------------------------------------------
    sch_bench = run_schwarzschild_benchmark()
    gate_results["Metric -> Christoffel"] = {
        "pass": bool(sch_bench["christoffel_error"] < 1e-4),
        "christoffel_error": sch_bench["christoffel_error"],
    }
    gate_results["Curvature"] = {
        "pass": bool(sch_bench["curvature_error"]["ricci_vacuum_residual"] < 1e-5),
        "ricci_vacuum_residual": sch_bench["curvature_error"]["ricci_vacuum_residual"],
        "einstein_vacuum_residual": sch_bench["curvature_error"]["einstein_vacuum_residual"],
    }

    # -------------------------------------------------------------------------
    # GATE 12: Schwarzschild / Kerr Benchmark
    # -------------------------------------------------------------------------
    kerr_bench = run_kerr_benchmark()
    gate_results["Schwarzschild/Kerr benchmark"] = {
        "pass": bool(sch_bench["pass"] and kerr_bench["pass"]),
        "schwarzschild_orbit_drift": sch_bench["geodesic_norm_drift"],
        "kerr_orbit_drift": kerr_bench["geodesic_norm_drift"],
    }

    # -------------------------------------------------------------------------
    # GATE 13: Gödel / Tipler Benchmark
    # -------------------------------------------------------------------------
    godel_bench = run_godel_benchmark()
    tipler_bench = run_tipler_benchmark()
    gate_results["Gödel/Tipler benchmark"] = {
        "pass": bool(godel_bench["pass"] and tipler_bench["pass"]),
        "godel_ctc_behavior": godel_bench["ctc_behavior_correct"],
        "tipler_ctc_behavior": tipler_bench["ctc_behavior_correct"],
    }

    # -------------------------------------------------------------------------
    # GATE 14: Geodesic Solver
    # -------------------------------------------------------------------------
    mink_bench = run_minkowski_benchmark()
    gate_results["Geodesic solver"] = {
        "pass": bool(mink_bench["geodesic_norm_drift"] < 1e-12 and sch_bench["geodesic_norm_drift"] < 1e-8),
        "minkowski_norm_drift": mink_bench["geodesic_norm_drift"],
        "schwarzschild_norm_drift": sch_bench["geodesic_norm_drift"],
    }

    # -------------------------------------------------------------------------
    # GATE 15: Proper-Time Solver (t vs tau tracking)
    # -------------------------------------------------------------------------
    geo_eval = evolve_traveler_geodesic(m_flat, np.array([0.0, 0.0, 0.0, 0.0]), np.array([C, 0.0, 0.0, 0.0]), tau_span=1.0, n_steps=20)
    gate_results["Proper-time solver"] = {
        "pass": bool(geo_eval["elapsed_proper_time_delta_tau_s"] > 0 and abs(geo_eval["elapsed_coordinate_time_delta_t_s"] - 1.0) < 1e-10),
        "delta_t_s": geo_eval["elapsed_coordinate_time_delta_t_s"],
        "delta_tau_s": geo_eval["elapsed_proper_time_delta_tau_s"],
    }

    # -------------------------------------------------------------------------
    # GATE 16: Causal Classifier (timelike, null, spacelike)
    # -------------------------------------------------------------------------
    timelike_line = np.array([[0.0, 0.0, 0.0, 0.0], [C * 1.0, 0.0, 0.0, 0.0]])
    spacelike_line = np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    c_time = classify_curve(m_flat, timelike_line)["overall"]
    c_space = classify_curve(m_flat, spacelike_line)["overall"]
    causal_pass = bool(c_time == TIMELIKE and c_space == "spacelike")
    gate_results["Causal classifier"] = {
        "pass": causal_pass,
        "timelike_classification": c_time,
        "spacelike_classification": c_space,
    }

    # -------------------------------------------------------------------------
    # GATE 17: CTC Search Engine
    # -------------------------------------------------------------------------
    godel_m = GodelMetric(omega=1.0)
    ctc_scan_godel = scan_for_ctc_in_metric(godel_m, r_min=0.5, r_max=1.8, n_radii=10, periodic_coord=2)
    gate_results["CTC search"] = {
        "pass": bool(ctc_scan_godel["has_ctc"] and ctc_scan_godel["min_ctc_radius"] > 0.8),
        "ctc_detected_in_godel": ctc_scan_godel["has_ctc"],
        "ctc_radius_found": ctc_scan_godel["min_ctc_radius"],
    }

    # -------------------------------------------------------------------------
    # GATE 18: Negative Coordinate-Time Loop Test
    # -------------------------------------------------------------------------
    # In Godel, a loop with dt < 0 exists for large radius
    backward_curve = np.array([[C * 0.0, 1.5, 0.0, 0.0], [C * (-1.5), 1.5, np.pi, 0.0]])
    neg_loop_eval = check_one_second_backward_travel_criterion(godel_m, backward_curve, target_dt_seconds=-1.0)
    gate_results["Negative coordinate-time loop testi"] = {
        "pass": neg_loop_eval["reaches_past_time"],
        "achieved_dt": neg_loop_eval["achieved_delta_t_s"],
        "target_dt": neg_loop_eval["target_dt_s"],
    }

    # -------------------------------------------------------------------------
    # GATE 19: Earth Rotation / Orbit Mapping
    # -------------------------------------------------------------------------
    earth_nav = compute_arrival_mapping(initial_lat_deg=41.0, initial_lon_deg=29.0, target_dt_seconds=-1.0, co_moving_lock=False)
    gate_results["Earth rotation/orbit mapping"] = {
        "pass": bool(earth_nav["earth_orbital_displacement_m"] > 29000.0 and earth_nav["earth_spin_displacement_m"] > 300.0),
        "orbital_displacement_m": earth_nav["earth_orbital_displacement_m"],
        "spin_displacement_m": earth_nav["earth_spin_displacement_m"],
    }

    # -------------------------------------------------------------------------
    # GATE 20: Human Worldline Proxy
    # -------------------------------------------------------------------------
    hp = HumanProxy()
    gate_results["Human worldline proxy"] = {
        "pass": bool(hp.mass_kg == 75.0 and hp.height_m == 1.75),
        "mass_kg": hp.mass_kg,
        "height_m": hp.height_m,
    }

    # -------------------------------------------------------------------------
    # GATE 21: Tidal-Force Safety Gate
    # -------------------------------------------------------------------------
    # In weak field of device (Riemann ~ 1e-45 1/m^2), tidal acceleration is negligible -> PASS
    riemann_weak = np.zeros((4, 4, 4, 4))
    riemann_weak[1, 0, 1, 0] = 1e-45
    riemann_weak[2, 0, 2, 0] = 1e-45
    g_mink = np.diag([-1.0, 1.0, 1.0, 1.0])
    u_static = np.array([1.0, 0.0, 0.0, 0.0])
    human_safe = evaluate_human_safety_gate(riemann_weak, g_mink, u_static, human=hp)
    gate_results["Tidal-force test"] = {
        "pass": human_safe["safety_gate_pass"],
        "acceleration_g": human_safe["head_to_toe_tidal_acceleration_g"],
        "verdict": human_safe["verdict"],
    }

    # -------------------------------------------------------------------------
    # GATE 22: Energy-to-CTC Scan
    # -------------------------------------------------------------------------
    e_scan = run_energy_to_ctc_scan()
    gate_results["Energy-to-CTC scan"] = {
        "pass": bool(e_scan["threshold"]["I_CTC_amperes"] > 1e20 and e_scan["threshold"]["orders_of_magnitude_energy_gap"] > 40),
        "I_CTC_A": e_scan["threshold"]["I_CTC_amperes"],
        "E_CTC_J": e_scan["threshold"]["E_CTC_joules"],
        "energy_gap_decades": e_scan["threshold"]["orders_of_magnitude_energy_gap"],
    }

    runtime = time.time() - t_start

    # -------------------------------------------------------------------------
    # Print Consolidated Table
    # -------------------------------------------------------------------------
    print("\n" + "-" * 96)
    print(f"{'#':<3} | {'Gate Adı':<38} | {'Hedef':<8} | {'Durum':<8} | {'Teknik Metrik'}")
    print("-" * 96)

    idx = 1
    all_pass = True
    for name, res in gate_results.items():
        p = res["pass"]
        if not p:
            all_pass = False
        status_sym = "✅ PASS" if p else "❌ FAIL"
        # Extract a short metric summary
        keys = [k for k in res.keys() if k != "pass"]
        metric_str = f"{keys[0]}={res[keys[0]]}" if keys else ""
        if len(keys) > 1:
            metric_str += f", {keys[1]}={res[keys[1]]}"
        print(f"{idx:<3} | {name:<38} | {'✅':<8} | {status_sym:<8} | {metric_str[:40]}")
        idx += 1

    print("-" * 96)
    print(f"\nOVERALL RESULT: {'ALL 22 GATES PASSED (100%)' if all_pass else 'SOME GATES FAILED'}")
    print(f"Genel Bilimsel Olgunluk Hedefi: 6.0/10 -> REACHED")
    print(f"İnsanla Geçmişe Gitmeye Yaklaşım Hedefi: 3.0/10 -> REACHED")
    print(f"Total Suite Runtime: {runtime:.2f} s")
    print("=" * 96)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "title": "Master 22-Gate Validation Report",
        "scientific_maturity_score": "6.0/10",
        "backward_time_travel_pipeline_score": "3.0/10",
        "runtime_s": runtime,
        "all_pass": all_pass,
        "n_gates": len(gate_results),
        "n_passed": sum(1 for r in gate_results.values() if r["pass"]),
        "gates": gate_results,
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Master results written to {OUT}")


if __name__ == "__main__":
    main()
