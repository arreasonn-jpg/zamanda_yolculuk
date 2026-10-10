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

## v2.0-alpha2: full dt_Gamma RHS

- `src/zt005/bssn_gamma.py` — conformal connection evolution
  (-2 Atilde^ij d_j alpha + 2 alpha [G^i_jk Atilde^jk
   - (2/3) gt^ij d_j K + 6 Atilde^ij d_j phi]).
- `tests/test_bssn_gamma.py` — 3 checks (Minkowski zero, finite, linear scale).

## v2.0-alpha3: gauge evolution

- `src/zt005/bssn_gauge.py` — 1+log lapse + Gamma-driver shift.
- `tests/test_bssn_gauge.py` — 5 checks.

## v2.0-alpha4: Z4c damping

- `src/zt005/bssn_z4c.py` — damping_dK, damping_dGamma, damping_energy.
- `tests/test_bssn_z4c.py` — 5 checks.

## v2.0-alpha5: unified driver

- `src/zt005/bssn_v2.py` — RK4 step on BSSNState with Gamma + gauge + Z4c.
- `scripts/long_run_v2.py` — 40-step run: growth sub-exponential,
  alpha_min = 1.000, all finite.
- `tests/test_bssn_v2.py` — 3 checks.

## v2.0-alpha5: unified driver

- `src/zt005/bssn_v2.py` — RK4 step on BSSNState with Gamma + gauge + Z4c.
- `scripts/long_run_v2.py` — 40-step run: growth sub-exponential,
  alpha_min = 1.000, all finite.
- `tests/test_bssn_v2.py` — 3 checks.

## v2.0-alpha5: real Hamiltonian + N=11 test

- `src/zt005/bssn_constraints.py` — grid-based H with true Ricci scalar.
- `tests/test_bssn_constraints.py` — 3 checks.
- N=7 vs N=11: growth per step 1.35 vs 1.58 — NOT converging.
- Root cause: missing terms in dAtilde (Gamma-dependent) and full Z4c.
- Next: v2.0-beta1 full dAtilde.

## v2.0-beta1: full physical Ricci + stable BSSN

- `src/zt005/bssn_ricci_full.py` — physical 3-Ricci scalar of
  g_ij = e^{4 phi} gtilde_ij, including phi-derivative terms.
- `src/zt005/bssn_full.py` — full RHS with covariant Hessian and D^2 alpha.
- `src/zt005/bssn_rk4_beta1.py` — RK4 driver with full RHS + real-H Z4c.

### Result (200-step run, N=7, dt=1e-4)

| Metric | Before (v2.0-alpha5) | After (v2.0-beta1) |
|---|---|---|
| H0 | 2.2e-12 (hayalet) | 6.08e-04 (fiziksel) |
| growth / 200 steps | ~1e+2 | **1.0009** |
| per-step growth | 1.585 | 1.0000045 |
| alpha_min | 1.0000 | 1.0000 |

Constraint is preserved to <0.1% over 200 RK4 steps. BSSN is now a
functionally stable 3+1 solver at skeleton level (shift=0).

## v2.0-beta2: shift evolution

- `src/zt005/bssn_shift.py` — Lie derivatives along beta for all fields.
- `src/zt005/bssn_rk4_beta2.py` — RK4 with shift evolution + gauge.
- `tests/test_bssn_shift.py` — 3 checks.
- `tests/test_bssn_rk4_beta2.py` — 3 checks.

### 100-step run (N=7, dt=5e-5)

| Metric | Value |
|---|---|
| H growth | **1.0001** |
| alpha_min | 1.0000 |
| beta_max | 2.3e-7 (linear growth) |
| all finite | True |

Shift evolves under Gamma-driver, constraint preserved.

## v2.0-beta3: Sommerfeld BC wired into RK4

- `src/zt005/bssn_bc.py` — sommerfeld_all for phi, K, alpha, gtilde,
  Atilde, Gamma, beta.
- `src/zt005/bssn_rk4_beta3.py` — RK4 driver with BC applied each step.
- `tests/test_bssn_bc.py` — 3 checks.
- `tests/test_bssn_rk4_beta3.py` — 3 checks.

### Note on integration test

BC is verified in unit tests. In the beta3 integration test the wave
does not reach the boundary within the run (0.1 time units vs 0.5
distance to boundary), so the BC is a no-op there. BC will trigger in
the beta4 long run (1000+ steps).
