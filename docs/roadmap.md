# Phase plan (locked)

Progression is gated: a phase starts only when the previous phase's
validation exists on record.

| Phase | Content | Status |
|---|---|---|
| 0 | Project infrastructure (repo layout, CI, docs, results policy) | ✅ this reorganization |
| 1 | Maxwell validation (MQS consistency) | ✅ ZT-005.x (full-wave deferred to ZT-007) |
| 2 | T_μν validation | ✅ EM-only T built; conservation diagnostic exists; device-source convergence RESOLVED via effective kernel smoothing a_wire_eff = 0.08 m (ZT-006.2.19/.20 FAIL → ZT-006.2.24) |
| 3 | Retarded linearized Einstein solver | ✅ solver implemented + manufactured-source validated; device-source convergence PASS at 5.8e-6 relative (ZT-006.2.25) — source chain numerically closed |
| 4 | Metric validation | ✅ benchmark metrics validated (known metrics) |
| 5 | Geodesic engine | ✅ RK4 integrator validated (norm, precession) |
| 6 | Known-CTC benchmark reproduction | ✅ Gödel + Tipler known answers reproduced; negative benchmarks clean |
| 7 | Device-generated spacetime search | ✅ completed — full chain I→(E,B)→T→h̄→g→CTC search for G1-G4; zero CTC candidates |
| 8 | Energy scaling (P_in/loss/stored/radiated/thermal, η, E_{1s}) | ✅ completed — energy_budget.py for all configs; P_in=84W, η~10⁻⁸-10⁻¹⁰, h~10⁻⁵⁰-10⁻⁵² |
| 9 | Earth + target-coordinate mapping (tensor level) | ✅ completed — Earth T^μν, linearized Schwarzschild h, combined metric; Earth dominates by 10⁴⁰-10⁴² |
| 10 | Engineering feasibility (coils, conductors, thermal, HV, power, backpack constraint optimizer) | ✅ completed — coil design (AWG8, R=80mm, N=10), thermal (forced air 50% duty), HV supply (PF correction, resonant cap), battery (500 Wh, 489 min), backpack budget; papers updated EN/TR, PDFs recompiled |
| 11 | Experimental apparatus (instrument-only; no human subjects) | 🔄 next |

## Immediate next moves (in order)

1. ~~**Conservation of the real device source.**~~ ✅ resolved via the
   ZT-006.2 chain (residuals tracked as regressions).
2. ~~**Retarded-solver convergence for the device source.**~~ ✅ ZT-006.2.24
   (instantaneous kernel, a_wire_eff) + ZT-006.2.25 (retarded solver
   PASS, 5.8e-6 relative).
3. ~~**Gauge audit of the device pipeline.**~~ ✅ completed — MQS non-conservation identified (scaled_residual ~ 0.85-1.13 vs 0.03 for conserved source).
4. ~~**Energy bookkeeping skeleton** (Phase 8 constants and budgets)~~ ✅ completed — energy_budget.py for all configs; P_in=84W, η~10⁻⁸-10⁻¹⁰, h~10⁻⁵⁰-10⁻⁵².
5. ~~**BSSN AMR frontier (v13/v14):**~~ ✅ v14.1 breakthrough with Brill-Lindquist; 50+ stable steps on L_c=200; moving-patch AMR working; projection reduces fine-grid band constraint by ~10-15%.
6. ~~**Phase 10: Engineering feasibility**~~ ✅ completed — coil design optimization, thermal analysis, HV/power supply, backpack constraint optimizer; papers updated EN/TR, PDFs recompiled.
7. **Phase 11: Experimental apparatus** — instrumentation design, safety protocols, data acquisition.

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
- ✅ P1a: Yakınsama hatası çözüldü (ZT-006.2.24, a_wire_eff=0.08 m);
  retarded-solver varyantı da PASS (ZT-006.2.25) → kaynak zinciri
  sayısal olarak kapandı.
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

- ✅ v2.0-beta5: Schwarzschild IC + 100-step vacuum evolution (H growth 1.0083)
- 🔷 v2.0-beta6: convergence study (N=9, 11, 13)

- ✅ v2.0-beta6: puncture IC + exterior convergence (H decreases 2.7x)
- 🔷 v2.0-beta7: Kerr benchmark (a>0)

- ✅ v2.0-beta7: slowly-rotating Kerr IC (a=0..0.5, growth 1.000)
- 🔷 v2.0-rc1: 1000-step long run + combined benchmark suite

- ✅ v2.0-rc1: 500-step + combined benchmark (growth <1.01)
- 🔷 v2.0: final release

## v2.0.x maintenance (in progress)

- ✅ physics_gap_summary.md updated with BSSN results
- ✅ README version table
- 🔷 paper/main.tex: BSSN section (draft, not submitted)
- 🔷 arXiv submission: ON HOLD (user decision)

## v2.0 scope note

v2.0 does NOT close the 43-order gap. It provides a stable non-linear
solver so future work can measure what happens near the CTC regime.

## v3.0 roadmap
- ✅ v3.0-alpha1: full Kerr KS Cartesian IC (|a| <= M)
- 🔷 v3.0-alpha2: moving-puncture formalism (gauge-IC consistent)
- 🔷 v3.0-alpha3: N-scaling with interpolated masks
- 🔷 v3.0-alpha4: horizon finding + mass extraction

- ✅ v3.0-alpha2: moving-puncture formalism (growth 1.0000)
- 🔷 v3.0-alpha3: Schwarzschild 500-step with moving puncture

- ✅ v3.0-alpha3: ADM mass module + moving-puncture long run
- 🔷 v3.0-alpha4: SC regular-core simulation (de Sitter + Israel matching)

- ✅ v3.0-alpha4: SC regular-core module (de Sitter + Israel junction)
- 🔷 v3.0-alpha5: SC core as BSSN IC + evolution test

- ✅ v3.0-alpha5: SC-like regular core BSSN IC (growth 1.0000)
- 🔷 v3.0-alpha6: SC core + Kerr (rotating SC)
- 🔷 v3.0-beta: full integration + release

- ✅ v3.0-alpha6: rotating SC (a=0..0.5, growth 1.0000)
- 🔷 v3.0-beta: consolidation + release

## v4.0 roadmap
- ✅ v4.0-alpha1: scalar field matter coupling (K ~ A^2, stable)
- 🔷 v4.0-alpha2: SC condensate as matter source
- 🔷 v4.0-alpha3: BSSN + SC core + matter (full integration)

- ✅ v4.0-alpha2: SC condensate as matter source (stable)
- 🔷 v4.0-alpha3: BSSN + SC metric + SC matter (full integration)

- ✅ v4.0-alpha3: BSSN + SC metric + SC matter (H growth 1.0000)
- 🔷 v4.0-beta: consolidation + release

## v5.0 roadmap
- ✅ v5.0-alpha1: PTUM field module (fuel + proxy)
- ✅ v5.0-alpha2: Tachyon condensation + NEC violation (rho_min = -15.6)
- 🔷 v5.0-alpha3: Long-run stability of negative-energy configuration

- ✅ v5.0-alpha3: diagnostics plot + long-run attempt (fails at ~100 steps)
- 🔷 v5.0-beta: consolidation (document 100-step stability limit)

- ✅ v5.0-alpha4: amplitude-dependent stability window identified
- 🔷 v5.0-beta: consolidation (document stability window)

## v6.0 roadmap
- ✅ v6.0-alpha1: complex scalar matter module
- ✅ v6.0-alpha2: RK4 on complex scalar + BSSN
- 🔷 v6.0-alpha3: constraint-preserving BC
- 🔷 v6.0-alpha4: AMR-lite (nested grid)

- ✅ v6.0-alpha3: two-level refinement operators
- 🔷 v6.0-alpha4: refined BSSN driver (coarse + fine grid coupled)

- ✅ v6.0-alpha4: refined scalar field (two-grid coupling)
- 🔷 v6.0-alpha5: refined BSSN driver (metric + matter on two grids)

- ✅ v6.0-alpha5: refined BSSN driver (two-grid coupling)
- 🔷 v6.0-alpha6: refined PTUM run (tachyon backreaction with refinement)

- ✅ v6.0-alpha6: refined PTUM 200-step stable (v5.0 crash resolved)
- 🔷 v6.0-beta: consolidation + release

- ✅ v6.0-alpha7: restrict mask-aware; PTUM discharge identified as transient
- 🔷 v6.0-beta: consolidate with transient-discharge interpretation

- ✅ v6.0-alpha8: Ford-Roman quantum inequality diagnostic
- 🔷 v6.0-alpha9: PTUM rho(t) sampled and tested against bound

- ✅ v6.0-alpha9: 1000-step refined run + Ford-Roman passed (no violation)
- 🔷 v6.0-beta: consolidate transient-discharge interpretation

- ✅ v6.0-alpha10: massive Ford-Roman (borderline; caveats documented)
- 🔷 v6.0-beta: consolidate; document open question on massive bound

## v6.1 roadmap
- ✅ v6.1-alpha1: unified framework document (SC + PTUM shared code)
- 🔷 v6.1-alpha2: speculative bridge note (tachyon mechanism common?)

- ✅ v6.1-alpha2: cross-scale phase transition bridge (speculative)
- 🔷 v6.1-beta: consolidation

## v6.2 roadmap
- ✅ v6.2-alpha1: Kerr + PTUM IC + single-grid baseline (near-crash)
- 🔷 v6.2-alpha2: refined Kerr + PTUM (2:1)

- ✅ v6.2-alpha2: refined Kerr + PTUM stable (100 steps, a=0..0.9)
- 🔷 v6.2-beta: consolidation + v6.2.0 release

- ✅ v6.2-alpha3: refined Kerr + PTUM 500-step stable (a=0.9)
- 🔷 v6.2-beta: consolidation + release

- ✅ v6.2-alpha4: Kerr + PTUM 1000-step stable (a=0.9)
- 🔷 v6.2-beta: consolidate (Kerr+PTUM as definitive PTUM-in-curved result)

- ✅ v6.3-alpha1: Kerr + SC regular core stable (50 steps, a=0..0.9)
- 🔷 v6.3-beta: consolidation + v6.3.0 release

- ✅ v6.3-alpha2: Kerr + SC + condensate full integration
- 🔷 v6.3-beta: consolidation + v6.3.0 release

- ✅ v6.3-alpha3: Kerr + SC 500-step (a=0.9, constraint improves)
- 🔷 v6.3-beta: consolidation + release

## v6.4 roadmap
- ✅ v6.4-alpha1: Kerr + SC + PTUM IC
- 🔷 v6.4-alpha2: full-stack evolution + 50-step test

## paper roadmap
- ✅ EN paper updated to v6.3
- ✅ TR scaffold + BSSN section
- 🔷 TR Introduction + Method (done)
- 🔷 TR Results, Discussion, Conclusion (pending)

- ✅ TR paper complete (all sections verified)
- 🔷 Final arXiv packaging (both EN and TR)

- ✅ paper v6.4 update (full-stack note EN+TR)
- ✅ arXiv v3 package assembled (982 KB)
- 🔷 Ready for arXiv upload (user decision pending)

- ✅ v6.4-alpha4: full stack SC + PTUM evolution (30-step stable)
- 🔷 v6.4-beta: consolidate + v6.4.0 release

- ✅ v6.5-alpha1: full stack 200-step long run (stable)
- 🔷 v6.5-alpha2: refined full stack (2:1 grid + SC + PTUM)

- ✅ v6.5-alpha2: refined full stack (100 steps, Kmax decreasing)
- 🔷 v6.5-beta: consolidate + v6.5.0 release

- ✅ v6.5-alpha3: refined full stack A=3.0 (NEC violation, 100-step stable)
- 🔷 v6.5-beta: consolidate + v6.5.0 release

- ✅ v6.6-alpha1: complex scalar + PTUM driver (activation needs seed)
- 🔷 v6.6-alpha2: psi seed + full activation test

- ✅ v6.6-alpha2: complex + PTUM equilibrium seed (NEC + stability)
- 🔷 v6.6-alpha3: refined complex + PTUM (A=3.0 with 2:1 grid)

- ✅ v6.6-alpha3: refined complex+PTUM (driver ok; equilibrium psi broken)
- 🔷 v6.6-beta: consolidate with noted limitation

- ✅ v7.0-alpha1: refined complex+PTUM NEC violation stable
- 🔷 v7.0-alpha2: consolidate as v7.0.0 release

- ✅ v7.1-alpha1: peak-preserving restrict (exploratory, two-way fails)
- 🔷 v7.1-alpha2: consolidate — one-way coupling remains canonical

- ✅ v7.1-alpha2: refined Kerr + complex + PTUM (NEC, 100-step)
- 🔷 v7.1-alpha3: Ford-Roman check on new config

- ✅ v7.1-alpha3: Ford-Roman violation detected at strong NEC
- 🔷 v7.1-beta: consolidate with QFT caveat

- ✅ v7.2-alpha1: transient discharge attempt (negative result documented)
- 🔷 v7.2-beta: consolidate as v7.2.0 release (negative result)

- ✅ v7.3-alpha1: QFT threshold scan (linear scaling, |rho|<1.27)
- 🔷 v7.3-alpha2: fine threshold scan + figure

- ✅ v7.3-alpha2: fine scan + threshold figure (ratio=0.795*|rho|)
- 🔷 v7.3-beta: paper integration + v7.3.0 release

- ✅ v7.4-alpha1: sub-QFT regime 4-regime structure identified
- 🔷 v7.4-alpha2: figure + paper integration

- ✅ v7.5-alpha1: discharge duration in QFT-OK window (frozen)
- 🔷 v7.5-beta: consolidate

- ✅ v7.5-alpha2: Kerr-invariance of QFT-OK window
- 🔷 v7.5-beta: consolidate

## v8.0 roadmap (new physics)
- ✅ v8.0-alpha1: Brans-Dicke matter module
- 🔷 v8.0-alpha2: BD evolution driver + GR limit test
- 🔷 v8.0-alpha3: Perihelion precession test (observational)
- 🔷 v8.0-alpha4: GW speed bound

- ✅ v8.0-alpha2: BD evolution driver + GR limit
- 🔷 v8.0-alpha3: perihelion precession (observational test)

- ✅ v8.0-alpha3: BD perihelion + GW speed observables
- 🔷 v8.0-alpha4: BD + PTUM combined (does BD affect QFT-OK window?)

- ✅ v8.0-alpha4: BD invariance of PTUM QFT window
- 🔷 v8.0-beta: v8.0.0 release + paper

## v8.1 roadmap (LQC)
- ✅ v8.1-alpha1: LQC Friedmann module
- 🔷 v8.1-alpha2: LQC-modified BSSN evolution (homogeneous)
- 🔷 v8.1-alpha3: bounce trajectory

- ✅ v8.1-alpha2: LQC scalar evolution
- 🔷 v8.1-alpha3: bounce trajectory (rho crosses rho_c)

- ✅ v8.1-alpha3: LQC bounce simulation (bounce detected)
- 🔷 v8.1-beta: paper integration + v8.1.0 release

## v8.2 roadmap
- ✅ v8.2-alpha1: Starobinsky f(R) scalar-tensor
- ✅ v8.2-alpha2: LQC + PTUM density comparison
- ✅ v8.2-alpha3: LQC QFT window invariance
- 🔷 v8.2-beta: v8.2.0 release

## v8.3 roadmap
- ✅ v8.3-alpha1: Hu-Sawicki f(R) dark energy
- ✅ v8.3-alpha2: Universal tensor-scalar QFT invariance
- 🔷 v8.3-beta: v8.3.0 release

## v9.0 roadmap
- ✅ v9.0-alpha1: Newman-Penrose GW extraction
- 🔷 v9.0-alpha2: GW extraction from BSSN evolution

## v9.2 roadmap
- ✅ v9.2-alpha1: GW extraction pipeline
- 🔷 v9.2-alpha2: Full BSSN run with GW extraction

- ✅ v9.2-alpha2: BSSN + GW extraction (null result for symmetric source)
- 🔷 v9.2-beta: document + v9.2.0 release

## v10.0 roadmap
- ❌ v10.0-alpha1/2: binary scalar GW (null result — physics)
- 🔷 v10.0-beta: document + release
- ✅ v10.1-alpha1/2: GW pipeline validated (strain extraction correct)
- 🔷 v10.1-beta: consolidate + paper integration
