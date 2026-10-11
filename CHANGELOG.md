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
