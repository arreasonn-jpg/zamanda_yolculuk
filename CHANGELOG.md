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
