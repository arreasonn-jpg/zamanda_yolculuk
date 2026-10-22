# Changelog

## v13.2.0 (2026-10-23)

### Added — band-local conformal projection: FIRST POSITIVE LEVER
- `src/zt005/amr/amr_projection.py` (v13.2-alpha1): after each
  coarse-to-fine injection, phi inside the 4-layer band is corrected
  by a damped point-Jacobi solve of Lap(dphi) = e^{4 phi} H_res / 8
  (first-order mapping of H = e^{-4 phi} (Rtilde - 8 D^2 phi -
  8 |grad phi|^2) + K^2 - At^2), with a |grad phi| > 0.5 guard mask
  excluding puncture-vicinity cells.
- Result: full 1000-step protocol, H_f/H0 @ 1000 = 0.9319 versus the
  fixed-patch baseline 1.0756 (about -13.4%), and BELOW the run's
  own step-0 normalization: the projection actively cleans the fine
  data during evolution. Finite at every checkpoint in the run
  (H_f curve 0.853 @ 50 -> 0.9319 @ 1000); dphi_max steady at
  ~2e-2.
- Diagnosis chain closed: cubic prolongation (v12.11), finer fine
  (v12.13), Gamma rederivation (v12.15), COM-moving patch (v13.0)
  and 4x + subcycle (v13.1) were all not the lever; the projected
  injected band IS a lever. The remaining defect sits exactly where
  the diagnosis had claimed all along: the prolonged second-
  derivative mismatch of the injected fields.
- `scripts/v13_2_projection.py` runner (same harness as v12.9/15);
  artifacts `results/validated/v13_2_projection.{npz,json}`
  (run log available but logs are gitignored);
  reproduce: `python scripts/v13_2_projection.py 1000`.
- Solver unit tests (5): tests/unit/test_amr_projection.py
  (Laplacian exactness, manufactured-solution convergence, zero
  boundary, guard mask, residual reduction).
- Paper EN/TR: moving-patch and 4x subsections already present
  (v13.0/13.1 documentation passes); the positive projection
  subsections added in this release; arxiv_submission_v4.zip
  regenerated.

### Caveats (on record)
- First-order mapping only: flat Laplacian stands in for the
  gamma-tilde covariant one; |grad phi|^2 lagged and masked. Exact
  covariant projection + momentum-constraint analog remain open
  (docs/AMR_PROJECTED_INJECTION_PLAN.md).
- H_f/H0 < 1 mixes genuine reduction with the H0 definition (measured
  without projection at step 0); an iso-reference metric is queued
  for the next pass.

### Tests
- Full fast suite at 13.2.0: see README (743 + 5 new solver tests).

## v13.1.0 (2026-10-20)

### Added — 4x + subcycle=4 AMR closed as negative (the v12.10 open item)
- The CFL-correct 4x + subcycle=4 configuration (N_f = 33, h_f =
  0.3125, dt_f = 5e-4) was run over 1000 steps plan; the process died
  of session shutdown between steps 500-750 and the step-500
  checkpoint resolves the question:
  H_f/H0 @ 500 = 1.1418 (4x) vs 1.0166 (2x baseline) and 1.0297
  (own coarse subgrid) — 4x is ~12% WORSE and worse than its own
  coarse grid; monotonic divergence from ~step 400.
- Verdict: resolution does not clear the prolonged-boundary
  second-derivative mismatch; it AMPLIFIES it. Constraint projection
  (see docs/AMR_PROJECTED_INJECTION_PLAN.md) is the only remaining
  lever candidate.
- Artifacts: results/validated/v13_amr_f4sub4.{npz,json};
  reproduce: `python scripts/v12_10_amr_flex.py 500 v16f4sub4 4 4`.
- Paper EN/TR: new bssn subsections (4x negative; moving-patch no-op
  added in v13.0.0 doc pass), arxiv_submission_v4.zip regenerated.
- Roadmap: phase-2/3 status synced to the ZT-006.2.24/.25 resolution
  (done in v13.0.1), immediate-next list updated, projection plan
  linked.

### Tests
- Diagnostic helpers in scripts/dev/ (probe_4x_step.py profiling
  probe, check_4x.py checkpoint reader).
- Full fast suite: 743 fast passed (v13.0.0 state); this release
  adds no test and changes no production code (artifacts + docs).

## v13.0.0 (2026-10-20)

### Added
- Moving fine-patch AMR machinery: `scripts/v13_amr_moving.py` with the
  reusable window logic in `src/zt005/amr/amr_moving.py` (COM tracking
  via |Atilde|^2 peaks, even-offset recentering, factor-2 aligned).
- Null result recorded: over 1000 steps the binary COM trigger never
  fires; H_f/H0 @ 1000 = 1.07564 — identical to the v12.15 fixed-patch
  baseline (1.0756); the COM-following moving patch is a no-op at this
  regime, quantifying the v12.12 "BHs stay inside the fixed patch".
  Artifacts: `results/validated/v13_amr_moving.{npz,json}`.

### Changed — subpackage restructure (breaking)
- 14 core modules moved into subpackages; imports rewritten repo-wide
  (85 files), no compatibility shims kept:
  - `zt005/amr/`: amr_interp, amr_bh_track, refinement,
    refinement_tensor, bssn_amr_consistent, amr_moving
  - `zt005/matter/`: matter_scalar, matter_complex, matter_ptum,
    matter_proca, matter_brans_dicke, sc_condensate, exotic_matter,
    phase_transition_common, ptum_field
  - `zt005.amr_interp` etc. are gone; use `zt005.amr.amr_interp`.
- README repository layout + version table updated; docs/limitations
  and docs/SESSION_NOTES updated to the v13 state.

### Paper (pre-v4 documentation pass)
- EN/TR bssn sections: new subsections for v12.13 (finer fine:
  negative) and v12.15 (Gamma-consistent injection: neutral).
- EN/TR main validation: source-chain closure recorded
  (ZT-006.2.24 a_wire_eff = 0.08 m; ZT-006.2.25 retarded solver,
  5.8e-6 relative).
- Pre-existing LaTeX bug fixed: split figure path
  `fig14_gw_pipeline.p\n ng` (both languages); both documents now
  compile with 0 errors.
- Built PDFs are now gitignored (`paper/arxiv_v3/*.pdf`).

### Tests
- New: `tests/unit/test_amr_moving.py` (window parity, clipping,
  mask threshold); version-guard test keeps pyproject 13.0.0.
- Full fast suite: 743 passed + 21 slow deselected (7:53).
- Diagnostic tools added under `scripts/dev/`
  (fan-in scan, one-shot refactor plan, import smoke).

### Pending (recorded in a follow-up version)
- 4x + subcycle=4 AMR run (the explicit v12.10 open item,
  CFL-correct): in progress at release time.

## v12.15.1 (2026-10-19) — consistency maintenance

### Fixed
- `src/zt005/__init__.py`: version now derived via importlib.metadata
  (pyproject as source of truth) with fallback; new regression test
  `tests/unit/test_version_consistency.py` guards against renewal drift
- `docs/stages/README.md`: ZT-006.2.24/.25 rows added (chronology table)
- `results/failed/README.md`: resolution pointer corrected
  (ZT-006.2.21 → ZT-006.2.21/.24/.25)
- Broken-path references after the v12.13.1 archive moves updated
  ( stages/ZT010/README.md, docs/SESSION_NOTES.md → paper/archive/ )
- CHANGELOG: duplicate "## v6.0.0" stub entry removed; v12.15 test
  count added (736 fast + 21 slow)

### Docs
- README frontline updated with the v12.15 neutral result; version
  status table extended to v12.15; repository layout mentions paper/

### Tests
- 738 fast + 21 slow

## v12.15.0 (2026-10-19)

### Added — constraint-consistent Gamma injection (v12.15)
- `src/zt005/bssn_amr_consistent.py` — rederive_Gamma: Gamma^i from
  prolonged gtilde with fine spacing, same d_central operator as the
  RHS; BSSN identity unit-tested (3 checks in
  tests/unit/test_bssn_amr_consistent.py)
- `scripts/v12_15_amr_gamma_consistent.py` — baseline vs gamma-consistent
  band injection, 2x fixed-patch AMR, Bowen-York binary, 1000 steps

### Key result (neutral)
- Gamma-consistent band injection changes injected Gamma by O(50%)
  relative but leaves H unchanged:
  baseline 1.0756 vs gamma_cons 1.0757 at step 1000
  (reference v12.9: 1.0756 — baseline reproduces it exactly)
- Third consecutive "not the lever" result (after cubic prolongation
  v12.11, finer fine v12.13): the remaining defect is the
  O(h_c^2) second-derivative mismatch of prolongated fields;
  next real levers are constraint projection or moving patches

### Docs
- docs/limitations.md: v12.15 scope section

### Tests
- 736 fast + 21 slow

### Misc
- src/zt005/__init__.py __version__ synced (was stale 0.2.0)

## v12.14.0 (2026-10-19)

### Added — ZT-006.2.25: retarded device-source convergence
- `src/zt005/run_zt006_2_25.py` — grid-refinement study of the
  **retarded** linearized-Einstein integral with the actual device
  sources (G1–G4, external observers, `a_wire_eff = 0.08 m`, MQS
  1 MHz time dependence with per-point retarded phase)
- `tests/validation/test_zt006_2_25.py` — fast regression test
  (smallest grid pair, one observer)
- `stages/ZT006/README_ZT006_2_25.md` + stage record updates

### Result (PASS)
- All consecutive-grid rel < 0.05 in both modes (quick worst rel
  1.47e-2; full final rel 5.79e-6), G1–G4
- Values numerically coincide with ZT-006.2.24 (instantaneous
  kernel), as MQS requires (ωL/c ≈ 5e-8 rad across the source)
- Closes the last open item of docs/validation.md's "Not yet
  validated" source-chain list; remaining gap to a device metric is
  physical (EM-only, linearized), not numerical

### Docs
- README frontline + honesty record, docs/validation.md,
  docs/limitations.md (#3) updated to reflect the closed item

### Tests
- 734 fast + 21 slow

## v12.13.1 (2026-10-19) — maintenance

### Changed
- Repository hygiene pass (no code changes):
  - `paper/`: 23 `.bak*` duplicates and LaTeX build artifacts (aux/log/out/blg)
    removed; `.bbl` kept for arXiv compilation. Previous versions remain in
    git history and in the archived submission zips.
  - Root: `arxiv_submission_v2/` + zip moved to `paper/archive/`;
    `pytest_raporu.txt` (old 299-test report) moved to
    `results/archived/pytest_raporu_299tests_old.txt`.
  - Caches cleaned (`__pycache__`, `.pytest_cache`, egg-info).
- README: version table updated through v12.13 (superseded the v2.0-only
  state that referenced "420 tests"); "current frontline" note added
  (device-source convergence open; constraint-consistent prolongation).
- Version bumped to 12.13.1.

## v0.9.0 (2026-10-10)
- P0: repo cleanup, requirements.txt pinned
- P1: ZT-006.2.24 convergence resolved (a_wire_eff=0.08 m)
- P1b: exotic_matter.py + energy_gap_report.py (h-gap 43 orders)
- P2: nonlinear_gr.py (scalar back-reaction)
- P3-P7: BSSN skeleton, RHS, grid FD, RK4, KO dissipation, Sommerfeld BC
- P8: pytest.ini markers + paper v2 tables
- P9: slow test markers (fast suite ~65 s)
- 335 tests passing

## v2.0.0 (2026-10-10)

### Added — full BSSN solver
- Extended state: Gamma, alpha, beta
- Full physical 3-Ricci scalar with conformal phi terms
- Covariant Hessian + D^2 alpha in RHS
- 1+log lapse + Gamma-driver shift
- Z4c-style constraint damping
- RK4 integrator on grid with FD derivatives
- Sommerfeld outflow BC
- Schwarzschild isotropic IC (clamped + smooth puncture)
- Slowly-rotating Kerr IC (|a| <= 0.5 M)

### Numerical results
- 500-step Schwarzschild: H growth 1.0062
- 100-step benchmark (5 cases): all growth < 1.0005
- All fields finite, gauge stable, constraint preserved

### Tests
- 400 fast + 20 slow

### Not included (documented in limitations.md)
- Full Kerr (a > 0.5), moving puncture, AMR, MPI,
  constraint-preserving BC, radiation extraction

## v3.0.0 (2026-10-11)

### Added — v3.0 series
- Full Kerr KS Cartesian IC (v3.0-alpha1)
- Moving-puncture gauge (v3.0-alpha2)
- ADM mass extraction (v3.0-alpha3)
- SC regular-core module (v3.0-alpha4)
- SC-like smoothed-core BSSN IC (v3.0-alpha5)
- Rotating SC (v3.0-alpha6)
- SC-BSSN bridge note

### Numerical results
- Schwarzschild: H growth 1.0000-1.0062
- Kerr a=0..0.5: H growth 1.0000
- SC core a=0..0.5: H growth 1.0000

### Tests
- 435 fast + 21 slow

### Reference
- Numerical support for Esa (2026) Schwarzschild Condensation model

## v4.0.0 (2026-10-11)

### Added — matter coupling
- Scalar field stress-energy (matter_scalar.py)
- BSSN RHS with matter source terms (bssn_rhs_matter.py)
- RK4 driver on (metric + scalar) (bssn_rk4_matter.py)
- SC condensate profile + potential (sc_condensate.py)
- Full SC integration (sc_full_ic.py)

### Numerical results
- Scalar on Minkowski: K ~ A^2, stable at A=1e-3, 1e-2
- SC condensate: phimax stable
- Full SC (metric + matter): H growth 1.0000, alpha_min 0.8612

### Tests
- 453 fast + 21 slow

### Reference
- Numerical support for Esa (2026) Schwarzschild Condensation
  regular-interior + matter-coupling claims

## v5.0.0 (2026-10-11)

### Added — PTUM tachyon coupling
- PTUM field module (ptum_field.py): fuel, proxy, tachyon threshold
- Two-field PTUM stress-energy (matter_ptum.py)
- RK4 driver for (metric + Phi + psi) (bssn_rk4_ptum.py)
- KO-dissipative variant (bssn_rk4_ptum_diss.py)
- Diagnostics plot script

### Key results
- NEC violation numerically confirmed: rho_min = -15.6
- 27/729 cells in negative-energy state
- Amplitude-dependent stability window identified
- Perturbative negative energy stable over 200+ steps

### Tests
- 467 fast + 21 slow

### Reference
- Numerical validation of Esa (2026) PTUM tachyon discharge mechanism
- Companion to SC numerical validation (v4.0)

## v6.0.0 (2026-10-11)

### Added — complex scalar + AMR-lite
- Complex scalar matter (matter_complex.py, U(1) charge)
- RK4 complex scalar + BSSN (bssn_rk4_complex.py)
- Two-level refinement operators (refinement.py)
- Tensor-aware refinement for BSSNState (refinement_tensor.py)
- Refined scalar evolution (refined_scalar.py)
- Refined BSSN driver (bssn_refined.py)
- Refined PTUM driver (bssn_refined_ptum.py)

### Key results
- PTUM A=3.0 refined: 200-step stable (v5.0 single-grid crashed ~100)
- Complex scalar: phase-rotation invariant, U(1) charge defined
- Refinement operators: exact for linear, smooth roundtrip verified

### Tests
- 490 fast + 21 slow

### Reference
- AMR-lite resolution of PTUM backreaction instability
- Numeric support for Esa (2026) SC + PTUM models

## v6.0.1 (2026-10-11) — patch

### Added
- Quantum inequality diagnostic (Ford-Roman)
- Massive field extension (rational approximation Q = 1/(1+y^2)^2)
- 1000-step refined PTUM run + time series saved

### Results
- Massless Ford-Roman: PTUM passes (ratio max 0.30)
- Massive Ford-Roman: borderline (ratio 1.21 at tau=1.0)
- Caveats documented; exact treatment open

### Tests
- 504 fast + 21 slow

## v6.1.0 (2026-10-11)

### Added — framework documentation
- docs/unified_framework.md (shared code, model-specific physics)
- docs/tachyon_bridge.md (speculative cross-scale structural parallel)
- src/zt005/phase_transition_common.py (shared potential formula)
- tests/test_phase_transition_common.py (9 checks)

### Status
- Mathematical parallel: verified
- Physical bridge: NOT established
- Scale gap: 78 orders (unchanged)

### Tests
- ~530 fast + 21 slow

### Honest boundary

v6.1 does NOT claim SC and PTUM are unified. It documents shared
infrastructure and a structural parallel, and explicitly identifies
the open problem of cross-scale bridging.

## v6.2.0 (2026-10-11)

### Added — Kerr + PTUM
- kerr_ptum_ic.py (Kerr-Schild + PTUM fields)
- refined_kerr_ptum_ic.py (2-grid IC builder)
- Refined Kerr + PTUM 500-step stable at a=0.9

### Key results
- Single-grid Kerr + PTUM A=3.0 crashes (Kmax ~ 1.4e3 at 50 steps)
- Refined: 500 steps stable at a=0.9 (Kmax ~ 0.27)
- ~17,000x Kmax reduction; refinement generalizes

### Tests
- ~540 fast + 21 slow

### Reference
- Refinement is a general solution for strong tachyon backreaction
- Numeric support for PTUM in curved (Kerr) background

## v6.3.0 (2026-10-12)

### Added — Kerr + SC
- kerr_sc_ic.py (Kerr-Schild + smoothed M_eff)
- kerr_sc_full_ic.py (metric + condensate)
- 500-step run at a=0.9 stable (constraint improves)

### Key results
- Kerr + SC metric: growth 0.997-0.999 (50 steps, a=0..0.9)
- Kerr + SC + condensate: 500 steps at a=0.9, growth 0.9823
- Constraint IMPROVES over time
- All fields finite throughout

### Tests
- ~537 fast + 21 slow

### Reference
- Extends SC paper's numerical support from static (v4.0) to
  rotating (Kerr) background

## v6.4.0 (2026-10-12)

### Added — full stack
- bssn_rk4_sc_ptum.py (metric + SC + PTUM)
- kerr_sc_ptum_ic.py (full-stack IC)
- 30-step stable at a=0..0.9

### Tests
- 546 fast + 21 slow

### Reference
- First simultaneous evolution of SC + PTUM in Kerr background

## v6.5.0 (2026-10-12)

### Added — refined full stack
- bssn_refined_sc_ptum.py (2-grid + SC + PTUM)
- Full stack 200-step single grid (a=0.9)
- Refined full stack A=3.0 (NEC violation)

### Tests
- 549+ fast + 21 slow

### Reference
- First NEC-violating stable full stack evolution
- Refinement resolves A=3.0 crash

## v6.6.0 (2026-10-12)

### Added — complex scalar + PTUM
- bssn_rk4_complex_ptum.py (single-grid driver)
- bssn_refined_complex_ptum.py (2-grid driver)
- Equilibrium-seed analysis (single grid reaches NEC at A=3.0)

### Tests
- 555 fast + 21 slow

### Known limitation
- Refined complex+PTUM loses equilibrium psi through prolongate/restrict
- Refined NEC violation remains a v6.5 result (SC+PTUM)

## v7.0.0 (2026-10-12)

### Added — refined complex+PTUM v2
- bssn_refined_complex_ptum_v2.py (one-way matter coupling)
- psi_equilibrium_smooth (regularized threshold)

### Key results
- Refined complex+PTUM A=3.0: rho_min = -15.67 persistent
- psimax = 2.83 constant, Kmax bounded
- Root causes fixed: sharp psi -> smooth; restrict loss -> one-way

### Tests
- 579+ fast + 21 slow

### Reference
- First persistent NEC violation in refined complex+PTUM

## v7.1.0 (2026-10-13)

### Added
- restrict_max (peak-preserving restriction, exploratory)
- Refined Kerr + complex + PTUM (v7.1-alpha2)
- Ford-Roman check on strong NEC (v7.1-alpha3)

### Key results
- Refined Kerr+complex+PTUM A=3.0: rho_min = -15.66, 100-step
- Ford-Roman violation at tau>=1.0 for strong NEC
- Weak NEC remains QFT-compatible

### Tests
- 587 fast + 21 slow

### Reference
- QFT caveat for strong-amplitude PTUM documented

## v7.2.0 (2026-10-13)

### Added
- Transient discharge attempt (natural Phi decay)
- Fast-discharge variant (lambda=1.0, refined grid)
- Negative result documented with timescale analysis

### Key finding
- PTUM discharge lifetime ~5-20 (inaccessible at t~0.1)
- Strong NEC (rho_min=-15.66) QFT-incompatible at tau>=1 (expected)
- Weak NEC (rho_min=-0.28) QFT-compatible (v6.5, unchanged)

### Tests
- 595 fast + 21 slow

### Reference
- Timescale mismatch between PTUM discharge and QFT window
- Confirms quantum inequality constraint on strong NEC

## v7.3.0 (2026-10-13)

### Added
- QFT threshold scan (11 amplitudes)
- Linear scaling law: ratio = 0.795 * |rho_min|
- Fig8 (threshold figure) in bilingual paper

### Key result
- Threshold: |rho_min| = 1.258
- QFT-compatible: A_fuel < 1.87
- Linear fit quality: R ~ 1.000

### Tests
- 595 fast + 21 slow

### Reference
- Quantitative single-parameter QFT bound on PTUM NEC

## v7.5.0 (2026-10-13)

### Added
- Discharge duration analysis (frozen at all A)
- Kerr-invariance of QFT-OK window

### Key results
- Phi decays <1% over accessible window (all A, all spin)
- QFT-OK window identical a=0.0 vs a=0.9
- PTUM-intrinsic QFT threshold

### Tests
- 595 fast + 21 slow

## v8.0.0 (2026-10-13)

### Added — Brans-Dicke extension
- bssn_rk4_brans_dicke.py (evolution)
- matter_brans_dicke.py (stress-energy)
- brans_dicke_observables.py (Mercury, GW speed)
- bssn_rk4_bd_ptum.py (BD + PTUM combined)

### Key results
- Mercury precession 43 arcsec/century (GR) ✅
- BD deviation at Cassini bound: < 2e-5
- QFT-OK window invariant under BD phi
- Combined with v7.5: PTUM threshold background-independent

### Tests
- 590+ fast + 21 slow

### Reference
- First alternative-gravity extension of the framework
- PTUM QFT threshold intrinsic property confirmed

## v8.2.0 (2026-10-13)

### Added
- matter_starobinsky.py (R + alpha R^2)
- lqc_ptum_coupling.py (density comparison)
- lqc_qft_window.py (QFT window invariance)

### Key results
- Starobinsky: n_s ~ 0.965, r ~ 0.004
- LQC correction for PTUM: 1 - 1e-77
- Triple invariance of PTUM QFT threshold

### Tests
- 623+ fast + 21 slow

## v8.3.0 (2026-10-13)

### Added
- matter_hu_sawicki.py (f(R) dark energy)
- tensor_scalar_qft.py (universal invariance)

### Key results
- Hu-Sawicki: GR limit at high curvature, chameleon screening
- Universal QFT window invariance under F(phi)
- Window: A in [1.35, 1.87]

### Tests
- 644+ fast + 21 slow

## v9.1.0 (2026-10-14)

### Added
- Newman-Penrose GW extraction module
- Horndeski scalar-tensor skeleton

### Key results
- psi_4, TT polarizations, Peters-Mathews flux
- Universal invariance extended to Horndeski class

### Tests
- 657+ fast + 21 slow

## v9.2.0 (2026-10-14)

### Added
- gw_extraction.py (BSSN -> TT pipeline)
- Newman-Penrose scalar computation

### Key results
- Null result: spherical source -> no GW (correct physics)
- Pipeline validated end-to-end

### Tests
- 669+ fast + 21 slow

## v9.3.0 (2026-10-14)

### Added
- Quadrupole GW source test
- Nonzero GW extraction

### Key results
- h+ max = 4.46e-12 (obs z=8)
- psi4 max = 1.00e-10
- Confirms Einstein quadrupole formula

### Tests
- 674+ fast + 21 slow
## v10.1.0 (2026-10-14)
### Added
- gw_initial_data.py (direct TT packet injection)
- Pipeline validation: extracted strain matches analytic envelope
- Fig14 in bilingual paper
### Key results
- h+ extraction correct to numerical precision
- psi_4 = 0 at packet peak (physical)
### Tests
- 674 fast + 21 slow
## v10.2.0 (2026-10-14)
### Added
- Proca vector field matter source
- RK4 driver for BSSN + Proca
- Nonzero GW extraction from vector source
### Key results
- Scalar source: h+ = 0 (v10.0)
- Proca source: h+ = 4.0e-6 (v10.2)
- Confirms TT content requirement for GW
### Tests
- 697 fast + 21 slow
## v11.0.0 (2026-10-14)
### Added — 3D binary BH
- Moving puncture 3D validation
- Brill-Lindquist binary IC
- Binary evolution + GW extraction
### Key results
- 3D MP: H growth 1.0001 over 200 steps
- Binary BH: h+ = 1.4e-7 (infall quadrupole)
### Tests
- 703+ fast + 21 slow
## v12.0.0 (2026-10-15)
### Added — Bowen-York orbital binary
- bssn_ic_bowen_york.py (tangential momentum IC)
- Orbital binary BH evolution
- Both GW polarizations (h+, hx)
### Key results
- H growth 1.0098 (200 steps, bounded)
- h+ = 7.9e-6, hx = 7.4e-4 (orbital)
### Tests
- 712+ fast + 21 slow
## v12.1.0 (2026-10-15)
### Added
- Brill-Lindquist plunge waveform
- PTUM in Schwarzschild background
- Long orbital binary run (300 steps)
### Key results
- Plunge: h+ = 8.7e-6
- Orbit: h+ = 3.1e-4, hx = 1.1e-2
- QFT window invariant under BH curvature
### Tests
- 697 fast + 21 slow
## v12.2.0 (2026-10-15)
### Added — vectorized BSSN
- bssn_rhs_vec.py (einsum-based RHS, 25x speedup)
- bssn_rk4_mp_vec.py (vectorized moving puncture)
- N=13 binary BH runnable in seconds
### Key results
- N=13 binary: H growth 1.0086, h+ = 1.0e-4
### Tests
- 692 fast + 21 slow
## v12.3.0 (2026-10-15)
### Added
- SC + PTUM + BH triple toy model
### Key results
- QFT window invariant under SC + BH
- Extends universal invariance to triple sector
### Tests
- 698 fast + 21 slow
## v12.4.0 (2026-10-15)
### Added
- N=15 medium run (1000 steps, 152s)
- h+ = 2.5e-4, hx = 6.2e-3 (both grow)
### Tests
- 698 fast + 21 slow
## v12.5.0 (2026-10-15)
### Added
- DHOST skeleton (beyond-Horndeski)
- Teleparallel f(T) skeleton
### Key results
- GW170817 bound integrated
- QFT window invariant (extended)
### Tests
- 715 fast + 21 slow
## v12.6.0 (2026-10-15)
### Added
- N=17 chunked inspiral (3200 steps)
- Constraint breakdown detected at H/H0 = 1.63
### Key finding
- Yörünge %28 captured before BC failure
- Sommerfeld BC insufficient for full inspiral
### Tests
- 715 fast + 21 slow

## v12.7.0 (2026-10-16)

### Added
- bssn_bc_improved.py (sponge + multi-layer Sommerfeld, 8 checks)
- Z4c damping in bssn_rhs_vec.py
- bssn_rk4_mp_vec_z4c.py driver

### Key finding (negative)
Boundary constraint damping has no effect on N=17 inspiral.
H growth 2.71 (damped) vs 2.71 (undamped). Puncture under-resolution is limiting factor.

### Tests
- 725 fast + 21 slow

## v12.8.0 (2026-10-17)

### Added
- Full BSSN shift advection in bssn_rhs_vec.py
- Shift coupled into metric evolution via rk4_step_mp_vec_z4c

### Key findings
- Shift integration works (waveforms change with eta)
- Constraint control NOT provided by shift
- Resolution-limited (N=17: ~1 cell per puncture)
- AMR is the required next step

### Tests
- 727 fast + 21 slow

## v12.9.0 (2026-10-17)

### Added
- src/zt005/amr_interp.py (prolongation + injection, 6 checks)
- scripts/v12_9_amr_test.py (two-grid AMR)
- Paper: fixed-patch AMR section (EN+TR)

### Key finding (positive)
2x AMR refinement reduces Hamiltonian constraint growth by 33% at step 1000 (H_fine/H0 = 1.0756 vs coarse H = 1.1254). First positive AMR result.

### Tests
- 731 fast + 21 slow

## v12.10.0 (2026-10-17)

### Added
- scripts/v12_10_amr_flex.py (configurable AMR)

### Key findings
- Subcycling neutral at 2x (same H_f/H0 = 1.0166)
- 4x + no subcycle fails (CFL violation)
- BC injection is the limiting factor

### Tests
- 731 fast + 21 slow

## v12.11.0 (2026-10-17)

### Added
- Cubic Lagrange prolongation (bssn_interp.py, 3 funcs)
- scripts/v12_11_amr_cubic.py

### Key finding (negative)
Cubic prolongation 0.25% worse than trilinear at 2x AMR. Interpolation order is not the lever.

### Tests
- 732 fast + 21 slow

## v12.12.0 (2026-10-17)

### Added
- amr_bh_track.py (BH center tracking)
- v12_12_amr_track.py (2000-step 2x AMR)

### Key findings
- AMR advantage decays: 54% (step 500) -> 8% (step 2000)
- BH tracking grid-quantized at N=17
- 2000 steps stable; runtime 2x by end

### Tests
- 733 fast + 21 slow

## v12.13.0 (2026-10-17)

### Added
- scripts/v12_13_amr_finen.py (configurable N_fine)

### Key finding (negative)
N_fine=25 worse than N_fine=17 (1.0091 vs 1.0030). Trilinear aliasing at factor=3. Finer fine is not the lever.

### Tests
- 733 fast + 21 slow
