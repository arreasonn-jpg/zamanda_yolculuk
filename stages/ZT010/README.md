# ZT-010 — Exotic-matter integration

**Status:** open, exploratory.

## Deliverables
- src/zt005/exotic_matter.py
- src/zt005/energy_gap_report.py
- tests/test_exotic_matter.py (7 checks)
- docs/engineering/energy_scaling.md
- results/validated/energy_gap_report.json

## Findings
| Model | log10(rho_ex/rho_EM) |
|---|---|
| Casimir (d=1um) | +26.6 |
| Morris-Thorne (b0=1m) | +38.8 |
| Tipler (w=1e6 rad/s) | +50.8 |

h-gap: 43 orders of magnitude.

Not a device. Numerical bound only.

## Nonlinear GR skeleton added

- `src/zt005/nonlinear_gr.py` — scalar source, second-order back-reaction,
  linear→nonlinear crossover.
- `tests/test_nonlinear_gr.py` — 5 checks.
- Crossover: h = 1 (order-unity), current h = 1e-46 → **46 orders**.
- Linear vs nonlinear ratio at h=1e-3: ~1e-6 (safe linear regime).
- Status: skeleton only; full BSSN solver not implemented.

## BSSN skeleton added

- `src/zt005/bssn.py` — ADM/BSSN containers, constraint residuals,
  Euler step, linear-regime check.
- `tests/test_bssn.py` — 6 checks.
- Status: skeleton only. No full BSSN evolution, no constraint damping.
- Purpose: lay out the 3+1 variable structure so future work can plug
  real RHS terms in `euler_step`.

## BSSN RHS added

- `src/zt005/bssn_rhs.py` — conformal Christoffel, Ricci, constraint
  damping (Z4c-inspired), RHS dict for `bssn.euler_step`.
- `tests/test_bssn_rhs.py` — 7 checks (flat limits + damping signs).
- Status: analytic pointwise RHS; grid + FD derivatives still external.

## BSSN evolution loop added

- `src/zt005/bssn_grid.py` — central FD derivatives on uniform 3D grid.
- `src/zt005/bssn_evolve.py` — Euler time loop + Hamiltonian constraint monitor.
- `tests/test_bssn_evolve.py` — 4 checks (FD linear, Minkowski stationary,
  constraint recorded).
- Status: pointwise Euler, no RK, no full dG (higher derivatives zero).

## RK4 + full dG added

- `src/zt005/bssn_rk4.py` — RK4 time integration.
- `bssn_grid.christoffel_from_grid` + `dG_from_grid` — full Gamma and
  its FD derivative fed into `ricci_tensor`.
- `bssn_evolve.evolve_bssn` now uses full dG (no more zero dummy).
- `tests/test_bssn_rk4.py` — 4 checks.
- Status: RK4 stable on Minkowski; constraint norm bounded.

## KO dissipation + Sommerfeld BC added

- `src/zt005/bssn_dissipation.py` — 6th-order Kreiss-Oliger + outflow BC.
- `tests/test_bssn_dissipation.py` — 4 checks.

## BSSN benchmark layer added

- `src/zt005/bssn_benchmark.py` — Schwarzschild isotropic initial data
  + conformal Ricci scalar via FD.
- `tests/test_bssn_benchmark.py` — 5 checks (flat gtilde, phi, alpha
  limits, Minkowski R=0, Schwarzschild R<1e-3).
- Fix: conformal Ricci scalar coefficient (Δφ + |∇φ|²), not (Δφ + ½|∇φ|²).

## Kerr + linear-regime cross-check added

- `src/zt005/bssn_kerr.py` — Kerr BL lapse + 3-metric, horizon radius.
- `src/zt005/bssn_linear_check.py` — hbar→metric, nonlinearity ratio.
- `tests/test_bssn_kerr.py` — 6 checks.

## Full BSSN RHS from ADM (partial) added

- `src/zt005/bssn_rhs_full.py` — dt_gtilde, dt_phi, dt_K, dt_Atilde,
  dt_Gamma, bssn_rhs_full (drop-in for euler_step).
- `tests/test_bssn_rhs_full.py` — 6 checks.
- Assumptions: alpha=1, beta=0; Atilde derivatives not yet wired in.

## 3D stability test added

- `src/zt005/bssn_stability.py` — Gaussian perturbation on Minkowski,
  evolve with RK4, constraint norm bounded.
- `tests/test_bssn_stability.py` — 4 checks.

## BSSN <-> retarded_hbar cross-check added

- `src/zt005/bssn_crosscheck.py` — weak-field h ~ 4 phi, consistency
  test against Newtonian 2M/r.
- `tests/test_bssn_crosscheck.py` — 4 checks.

## Grid-aware full BSSN RHS (P-A1)

- `src/zt005/bssn_rhs_grid.py` — full dphi/dgt/dK/dAtilde with FD
  derivatives of phi (grad + Hessian) and Atilde wired in.
- `tests/test_bssn_rhs_grid.py` — 4 checks.

## Grid-aware RK4 driver (P-A2)

- `src/zt005/bssn_rk4_grid.py` — RK4 with grid RHS (Hessian phi,
  grad Atilde, full Ricci).
- `tests/test_bssn_rk4_grid.py` — 3 checks.

## v2.0-alpha1: extended BSSN state

- `src/zt005/bssn_state.py` — `BSSNState` dataclass with phi, gtilde, K,
  Atilde, Gamma, alpha, beta; `minkowski_state(N)` helper.
- `tests/test_bssn_state.py` — 4 checks.
