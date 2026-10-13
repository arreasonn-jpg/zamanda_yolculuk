# Changelog

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

## v6.0.0 (2026-10-11)

### Added — complex scalar + AMR-lite
(see previous entry)

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
