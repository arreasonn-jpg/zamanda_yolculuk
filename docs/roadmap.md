# Phase plan (locked)

Progression is gated: a phase starts only when the previous phase's
validation exists on record.

| Phase | Content | Status |
|---|---|---|
| 0 | Project infrastructure (repo layout, CI, docs, results policy) | ✅ this reorganization |
| 1 | Maxwell validation (MQS consistency) | ✅ ZT-005.x (full-wave deferred to ZT-007) |
| 2 | T_μν validation | 🔶 partial: EM-only T built; conservation diagnostic now exists but device-source convergence FAILED (ZT-006.2.19/.20) |
| 3 | Retarded linearized Einstein solver | 🔶 solver implemented + manufactured-source validated; device-source convergence study pending |
| 4 | Metric validation | ✅ benchmark metrics validated (known metrics) |
| 5 | Geodesic engine | ✅ RK4 integrator validated (norm, precession) |
| 6 | Known-CTC benchmark reproduction | ✅ Gödel + Tipler known answers reproduced; negative benchmarks clean |
| 7 | Device-generated spacetime search | ❌ not started — gated on Phases 2/3 convergence |
| 8 | Energy scaling (P_in/loss/stored/radiated/thermal, η, E_{1s}) | ❌ not started |
| 9 | Earth + target-coordinate mapping (tensor level) | ❌ kinematic placeholder only |
| 10 | Engineering feasibility (coils, conductors, thermal, HV, power, backpack constraint optimizer) | ❌ not started |
| 11 | Experimental apparatus (instrument-only; no human subjects) | ❌ not started |

## Immediate next moves (in order)

1. **Conservation of the real device source.** Evaluate ∂_μ T^{μν} for
   the MQS device fields (including the wire currents explicitly) and
   make the residual a tracked regression.
2. **Retarded-solver convergence for the device source.** The failure
   mode of ZT-006.2 must not repeat: grid-refinement study with
   analytic/semi-analytic sub-checks before any metric claim.
3. **Gauge audit of the device pipeline.** ∂_μ h̄^{μν} for the device
   retarded solution.
4. **Energy bookkeeping skeleton** (Phase 8 constants and budgets), so
   that the h ~ 10⁻⁴⁶ scale can be translated into an explicit energy
   requirement.

## Explicitly deferred

* C/C++ kernels and GPU acceleration — not before the Python baseline is
  right (a fast wrong answer is worse than a slow right one).
* Docker / packaging polish — after Phase 3 closure.
* Full-wave Maxwell (ZT-007) — when frequency scans demand it.

## Human-experiment policy

Any human experiment is out of scope for all listed phases and remains
so until an independent ethics/safety review concludes otherwise.

## Update — P0/P1 tamamlandı

- ✅ P0: Kök dizin temizlendi, `requirements.txt` üretildi, 305 test geçiyor.
- ✅ P1a: Yakınsama hatası çözüldü (ZT-006.2.24, a_wire_eff=0.08 m).
- ✅ P1b: `exotic_matter.py` + `energy_gap_report.py` eklendi.
- 🔷 P2: Doğrusal olmayan rejim (BSSN) henüz yok.

- ✅ P3: `bssn.py` (3+1 skeleton) + `nonlinear_gr.py` (scalar back-reaction).
- 🔷 P4: Full BSSN RHS + constraint damping (still open).

- ✅ P4: BSSN RHS pointwise (Christoffel, Ricci, constraint damping).
- 🔷 P5: Grid + finite-difference derivatives + constraint monitoring loop.

- ✅ P5: Grid FD + evolution loop + constraint monitor.
- 🔷 P6: RK4 time integration + full dG (Christoffel derivatives).

- ✅ P6: RK4 + full dG (Christoffel derivatives into Ricci).
- 🔷 P7: Kreiss-Oliger dissipation + boundary conditions.

- ✅ P7: Kreiss-Oliger dissipation + Sommerfeld boundary.
- 🔷 P8: pytest markers for slow tests.

- ✅ P-C1: BSSN benchmark (Schwarzschild isotropic, Ricci ~ 0).
- 🔷 P-C2: Kerr benchmark + linear-regime cross-check with ZT-006.3.

- ✅ P-C2: Kerr + linear-regime cross-check.
- 🔷 P-C3: BSSN ↔ retarded_hbar numeric cross-check.

- ✅ P10: Full BSSN RHS from ADM (lapse/shift fixed).
- 🔷 P11: 3D stability + constraint-preserving BC.

- ✅ P11: 3D stability (small Gaussian, RK4, bounded constraint).
- 🔷 P12: BSSN <-> retarded_hbar numeric cross-check.

- ✅ P12: BSSN <-> retarded_hbar cross-check.
- 🔷 P13: final "gap to Terminator-5" physics doc + v1.0 tag.

- ✅ P-A1: grid-aware full BSSN RHS (Hessian phi + grad Atilde).
- 🔷 P-A2: RK4 driver with grid RHS + long-run stability.

- ✅ P-A2: RK4 grid driver.
- 🔷 P-A3: long-run stability (100+ steps).

## v2.0 roadmap (Tam BSSN solver)
- ✅ v2.0-alpha1: extended state container (Gamma, alpha, beta)
- 🔷 v2.0-alpha2: full dt_Gamma RHS
- 🔷 v2.0-alpha3: 1+log lapse + Gamma-driver shift
- 🔷 v2.0-alpha4: Z4c damping
- 🔷 v2.0-alpha5: long-run stability (1000+ steps)
- 🔷 v2.0-beta: Schwarzschild + Kerr benchmarks
- 🔷 v2.0: release

- ✅ v2.0-alpha2: full dt_Gamma RHS
- 🔷 v2.0-alpha3: 1+log lapse + Gamma-driver shift

- ✅ v2.0-alpha3: 1+log lapse + Gamma-driver shift
- 🔷 v2.0-alpha4: Z4c constraint damping

- ✅ v2.0-alpha4: Z4c damping layer
- 🔷 v2.0-alpha5: unified v2.0 driver + long-run stability

- ✅ v2.0-alpha5: unified v2.0 driver + 40-step run
- 🔷 v2.0-alpha6: Sommerfeld outflow BC in v2.0 driver

- ✅ v2.0-alpha5: unified v2.0 driver + 40-step run
- 🔷 v2.0-alpha6: Sommerfeld outflow BC in v2.0 driver

- ✅ v2.0-alpha5 + realH constraint grid
- 🔷 v2.0-beta1: full dAtilde with Gamma terms

- ✅ v2.0-beta1: full physical Ricci -> stable BSSN (200-step H growth 1.0009)
- 🔷 v2.0-beta2: shift vector evolution (beta != 0)
- 🔷 v2.0-beta3: Sommerfeld BC in full driver
- 🔷 v2.0-beta4: long-run 1000+ steps

- ✅ v2.0-beta2: shift vector evolution (H growth 1.0001)
- 🔷 v2.0-beta3: Sommerfeld BC in full driver

- ✅ v2.0-beta3: Sommerfeld BC in full driver (unit-verified)
- 🔷 v2.0-beta4: long-run 1000+ steps (wave reaches boundary)

- ✅ v2.0-beta4: BC-engaged run (edge constant, BC works)
- 🔷 v2.0-beta5: Schwarzschild isotropic constraint-satisfying IC
