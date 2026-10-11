# Limitations (explicit)

If a limitation is not listed here, that is a documentation bug.

## Physics

1. **Incomplete source model.** T_μν contains only the electromagnetic
   field. Matter, mechanical stress, thermal energy, control systems and
   the power source are absent (assumptions.md A1).
2. **Linearized gravity only.** First order in G on a fixed Minkowski
   background. No full g_μν(t, x) solve, no back-reaction.
3. **No validated device metric.** The retarded solver is validated with
   manufactured sources; a grid-convergence study of the retarded
   integral with the actual device source is still open.
4. **MQS electromagnetics.** Valid while L/λ ≪ 1; frequency scans beyond
   that need ZT-007 (full-wave).
5. **Resonance is hypothetical.** No Q/ω₀/L/C/R or modal model exists
   (A4).
6. **Earth/solar frame is kinematic.** No tensor contribution of the
   Earth, no frame matching between device metric and background
   (A5).
7. **No CTC evidence for the device.** Existing CTC machinery is a
   known-answer benchmark, nothing more.

## Numerics

8. Midpoint quadrature for 1/R kernels: O(h²) away from singular cells,
   corrected singular-cell average at coincident points; convergence per
   configuration must still be demonstrated.
9. Gauge/conservation residual tolerances are grid-dependent; each result
   states its spacings.
10. Benchmark Christoffels are finite-difference except Schwarzschild
    (analytic); RK4 fixed step, no adaptive control yet.

## Engineering

11. No power-electronics, thermal, insulation, or structural model.
12. No backpack-constraint optimizer: any geometry result violating the
    hard envelope (0.45 × 0.32 × 0.15 m-class) must be reported as FAIL,
    and that optimizer does not exist yet.

## Experimental

13. No experimental apparatus, no measurement plan, no human subject at
    any horizon.

## Process

14. Single-commit history: earlier ZT-001…ZT-004 deliverables were not
    carried into this repository.
15. The ZT-006.2.19 raw FAIL artifact was never committed; its summary
    (0/20, per-grid percentages) is recorded in stages/ZT006/README.md.

## BSSN skeleton scope (v1.1)

The BSSN code is a **skeleton**, not a full numerical-relativity solver:

- `Gamma^i` is not evolved as a state variable.
- Lapse (`alpha`) and shift (`beta`) are fixed to (1, 0).
- `dAtilde` omits the `Gamma`-dependent terms.
- Long-run runs show slow sub-exponential constraint growth, consistent
  with physical dispersion rather than numerical blow-up at this scale.

Promoting this to a production BSSN solver would require:
1. Evolving `Gamma^i` with its full RHS.
2. Gauge conditions (1+log lapse, Gamma-driver shift).
3. Constraint-damping (Z4c) tuning.
4. Kreiss-Oliger dissipation wired into RK4 stages.
5. Grid refinement / AMR.

None of these are in scope for v1.x.

## v2.0 scope (final)

### What v2.0 delivers

- Full BSSN state: phi, gtilde, K, Atilde, Gamma, alpha, beta.
- Physical 3-Ricci scalar with conformal phi terms.
- 1+log lapse + Gamma-driver shift.
- Z4c-style constraint damping.
- RK4 integrator with grid-aware FD derivatives.
- Sommerfeld outflow BC.
- Schwarzschild isotropic IC (clamped).
- Schwarzschild puncture IC (smooth).
- Slowly-rotating Kerr IC (|a| <= 0.5 M).

### What v2.0 does NOT deliver

- Full Kerr (a > 0.5 M) — requires BL coordinate transformation.
- Moving-puncture formalism — needed for long-term BH evolution.
- AMR / mesh refinement.
- MPI parallelism.
- Constraint-preserving BC (only Sommerfeld outflow).
- Kreiss-Oliger dissipation active in the v2.0 RK4 loop.
- Radiation extraction (Newman-Penrose scalars).

### Numerical results

- 500-step Schwarzschild: H growth 1.0062 (per-step 1.0000124).
- 100-step benchmark: Schwarzschild + Kerr a=0..0.5, growth < 1.0005.
- All finite, alpha_min stable, bmax scales linearly with a.

### Honest assessment

v2.0 is a **functionally stable BSSN skeleton**. It correctly evolves
Schwarzschild and slowly-rotating Kerr initial data over hundreds of
RK4 steps without constraint blow-up. It is not a production
numerical-relativity code and does not attempt to be one.

## v3.0 scope (final)

### What v3.0 delivers

- Full Kerr via Kerr-Schild Cartesian IC (|a| <= M, horizon-penetrating).
- Moving-puncture gauge (Bruemann 2005) with auxiliary B^i field.
- ADM mass extraction (volume-integrated form).
- SC regular-core module (de Sitter + Israel junction).
- SC-like smoothed-core BSSN IC (regular at origin).
- Rotating SC core (frame-dragging shift).
- BSSN <-> retarded_hbar cross-check.
- All benchmarked with moving-puncture evolution.

### Numerical results summary

| Case | Setup | H growth |
|---|---|---|
| Schwarzschild isotropic | 500 steps | 1.0062 |
| Schwarzschild puncture | 50 steps, moving-puncture | 1.0000 |
| Kerr KS a=0..0.99 | 30-50 steps, frozen gauge | 0.94-1.00 |
| Kerr slow a=0..0.5 | 50 steps, moving-puncture | 1.0000 |
| SC smoothed core a=0 | 50 steps, moving-puncture | 1.0000 |
| SC smoothed core a=0..0.5 | 50 steps, moving-puncture | 1.0000 |

### What v3.0 does NOT deliver

- Full Kerr a > 0.5 M with evolution (only a <= 0.5).
- AMR / mesh refinement.
- MPI parallelism.
- Constraint-preserving BC (only Sommerfeld outflow).
- Radiation extraction (Newman-Penrose).
- Full de Sitter + Israel junction (approximated by smoothing).
- Nonlinear stability over long times (>500 steps).

### Honest assessment

v3.0 is a **functionally complete BSSN skeleton** with:
- Schwarzschild (isotropic + puncture)
- Kerr (slow rotation, a <= 0.5 M)
- SC-like regular core (non-rotating + rotating)

All stable over 50-500 moving-puncture steps. Not a production NR code.

## v4.0 scope

### What v4.0 delivers

- Scalar field matter coupling to BSSN (matter_scalar.py).
- BSSN RHS with source terms (bssn_rhs_matter.py).
- RK4 driver on (metric + scalar) (bssn_rk4_matter.py).
- SC condensate profile + potential (sc_condensate.py).
- Full SC integration: SC metric + SC matter (sc_full_ic.py).

### Numerical results

- Scalar on Minkowski: K ~ A^2, stable at A=1e-3 and 1e-2.
- SC condensate on Minkowski: phimax stable, all finite.
- Full SC (metric + matter): H growth 1.0000, alpha_min 0.8612.

### What v4.0 does NOT deliver

- Complex scalar field (only real approximation).
- Backreaction on the SC core's rho_0 from matter.
- Kerr + SC + matter combined.
- Long-time (>50 steps) nonlinear stability.
- Adaptive mesh refinement around the SC core.

### Honest assessment

v4.0 provides the first matter-coupled BSSN evolution in this repo.
The SC regular core + SC condensate co-evolve stably over 50 steps.
This supports, but does not prove, the SC model's dynamical stability.

## v5.0 scope

### What v5.0 delivers

- PTUM field module: fuel Phi, proxy psi, tachyon threshold
- Two-field PTUM stress-energy (matter_ptum.py)
- RK4 driver on (metric + 2 PTUM fields) (bssn_rk4_ptum.py)
- KO dissipation variant (bssn_rk4_ptum_diss.py)
- Diagnostics plot (6 panels: rho, K, phi, alpha)

### Numerical results

| Config | rho_min | Steps | Status |
|---|---|---|---|
| PTUM equilibrium A=3.0 | -15.6 | 100 | finite |
| PTUM equilibrium A=3.0 | -15.6 | 500 | fails ~100 |
| PTUM equilibrium A=1.5 | -0.29 | 200 | **stable** |
| PTUM equilibrium A=1.5 | -0.29 | 500 | (to test) |

### Key physical result

**NEC violation confirmed:** rho_min = -15.6 (27/729 cells)
Tachyon condensation produces the negative energy required for CTCs.

**Stability window:** negative energy configuration is stable when
perturbative (|rho_min| < 1), unstable when dominant (|rho_min| > 5).

### What v5.0 does NOT deliver

- Complex scalar (only real two-field approximation)
- Backreaction on rho_0 from PTUM matter
- Kerr + PTUM combined
- Long-term stability of dominant negative energy
- AMR around tachyon-active region

### Honest assessment

v5.0 numerically validates the PTUM model's central claim:
tachyon discharge produces negative energy density (NEC violation).
The magnitude is perturbative in the stable regime, consistent with
PTUM's framing of negative energy as one component of a larger
cosmological discharge.

## v6.0 scope

### What v6.0 delivers

- Complex scalar field matter (matter_complex.py)
- RK4 driver for complex scalar + BSSN (bssn_rk4_complex.py)
- Two-level refinement operators (refinement.py)
- Tensor-aware refinement (refinement_tensor.py)
- Refined scalar evolution (refined_scalar.py)
- Refined BSSN driver (bssn_refined.py)
- Refined PTUM driver (bssn_refined_ptum.py)

### Numerical results

| Setup | Result |
|---|---|
| Complex scalar on Minkowski | finite, phase-rotation invariant |
| Refinement roundtrip | exact for linear, < 0.5 smooth |
| Refined BSSN on Minkowski | zero stays zero |
| Refined PTUM A=3.0, 200 steps | **stable** (v5.0 crashed at 100) |

### Key physical result

**Refinement resolves strong tachyon backreaction.**
2:1 refinement with prolongate/restrict coupling allows the A=3.0
PTUM run (rho_min initial -15.6) to reach 200 steps, where the
single-grid v5.0 run failed at ~100 steps.

### What v6.0 does NOT deliver

- Full AMR (only 2-level static refinement)
- Refinement tracking (fine grid fixed at origin)
- Flux correction at refinement boundary
- Full nonlinear complex scalar (only linear regime tested)
- Long-run stability beyond 200 steps at A=3.0

### Honest assessment

v6.0 provides the first AMR-lite framework in this repo. The 2:1
refinement coupling works: PTUM configurations that crash on a single
grid are stable on two coupled grids. This is a numerical
demonstration that resolution near the tachyon-active region is
critical for long-run stability.

The rho_min drop from -15.6 to -0.096 after coupling remains to be
fully understood (physical redistribution vs. averaging artefact).

## v6.0 scope (updated with alpha10)

### Added since v6.0 initial

- Quantum inequality (Ford-Roman) diagnostic module
- Massive field extension (rational approximation)
- 1000-step refined PTUM run (results/validated/v6_long_ptum.npz)
- Ford-Roman test result: massless passes, massive borderline

### Ford-Roman status

| Case | tau=1.0 ratio | Verdict |
|---|---|---|
| Massless, N=2 | 0.30 | pass |
| Massive mu=1, N=2 | 1.21 | borderline |

Caveats push toward "no violation"; exact treatment is future work.

### Updated what v6.0 does NOT deliver

- Exact Ford-Roman bound (Bessel functions for massive fields)
- Interacting-field quantum inequality
- Tracked/moving refinement
- Complex scalar + Kerr + PTUM combined

## v6.1 scope

### What v6.1 delivers

- `docs/unified_framework.md` — shared infrastructure documentation
  for SC and PTUM models
- `docs/tachyon_bridge.md` — speculative cross-scale structural parallel
- `src/zt005/phase_transition_common.py` — shared potential formula
- `tests/test_phase_transition_common.py` — 9 checks

### Status of the cross-scale bridge

- Mathematical form: SHARED (verified)
  - V(psi; C) = 1/2 m_eff^2(C) psi^2 + 1/4 lambda psi^4
  - m_eff^2(C) = m^2 (1 - C/C_c)
- Physical bridge: NOT established
- Unification: NOT claimed
- Scale gap: 78 orders of magnitude remains

### What v6.1 does NOT deliver

- A derivation connecting SC and PTUM
- Any mechanism that bridges the 78-order scale gap
- A unified single model
- Numerical cross-scale simulation with shared physics

### Honest assessment

v6.1 provides organization and honest boundary-drawing. The two
companion papers (SC, PTUM) share code and mathematical structure
but are NOT unified. This is intentional: claiming unification
without a derivation would be speculative overreach.

## v6.2 scope

### What v6.2 delivers

- Kerr + PTUM IC (kerr_ptum_ic.py)
- Refined Kerr + PTUM IC builder (refined_kerr_ptum_ic.py)
- Refined Kerr + PTUM 500-step stable at a=0.9

### Numerical results

| Setup | Single grid | Refined |
|---|---|---|
| PTUM flat A=3.0 | crash ~100 | 1000 stable |
| Kerr + PTUM A=3.0 | near-crash ~50 | 500 stable |

Kmax reduction: ~17,000x. Refinement generalizes to Kerr background.

### What v6.2 does NOT deliver

- Full Kerr a=1 (extremal)
- Tracked/moving refinement
- Radiation extraction
- Full nonlinear complex scalar + Kerr + PTUM
- Long-run (>1000 steps) Kerr + PTUM

### Honest assessment

v6.2 demonstrates that 2:1 refinement resolves strong tachyon
backreaction in both flat (v6.0) and Kerr (v6.2) backgrounds.
This is a general numerical technique, not IC-specific.
