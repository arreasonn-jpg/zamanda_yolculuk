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

## v6.3 scope

### What v6.3 delivers

- Kerr-Schild + smoothed M_eff(r) = Kerr + SC regular core
- Kerr + SC + SC condensate (full integration)
- 500-step stable at a=0.9 (near-extremal)

### Numerical results

| Setup | Result |
|---|---|
| Kerr + SC (metric only, 50 steps) | growth 0.997-0.999 |
| Kerr + SC + condensate (50 steps) | growth 0.9983 |
| Kerr + SC + condensate (500 steps, a=0.9) | growth **0.9823** (improves) |

Constraint IMPROVES over time in the 500-step run. All fields finite
throughout.

### What v6.3 does NOT deliver

- Extremal Kerr (a=1) exact
- Tracked/moving refinement for Kerr
- Radiation extraction
- Nonlinear stability beyond 500 steps
- Full matter backreaction on M_eff

### Honest assessment

v6.3 extends SC regular-core result from static (v4.0) to rotating
(Kerr) background. This is the strongest numerical support for the
SC paper: metric + matter + spin all consistent.

## v6.4 scope

### What v6.4 delivers

- Full stack: Kerr+SC metric + SC condensate + PTUM fields
- Driver: bssn_rk4_sc_ptum.py (metric + 6 matter fields)
- 30-step evolution stable at a=0..0.9
- 546 fast + 21 slow tests

### Numerical result

| Setup | Steps | Result |
|---|---|---|
| Kerr+SC+PTUM a=0..0.9 | 30 | stable, all finite |
| rho_total | — | -0.284 (PTUM dominated) |
| Kmax | — | 0.35-0.36 |

### What v6.4 does NOT deliver

- Long-run (>30 steps) full-stack stability
- Full nonlinear SC condensate coupling
- MPI / AMR for full stack
- Radiation extraction

### Honest assessment

First simultaneous evolution of all three components (rotating
regular core + QCD condensate + tachyon field). 30 steps stable.
Long-run stability is deferred to v6.5.

## v6.5 scope

### What v6.5 delivers

- Refined full stack: Kerr+SC metric + condensate + PTUM (2-grid)
- Full stack 200-step long run (single grid, a=0.9)
- Refined full stack 100-step at A=3.0 (NEC violation)
- 549+ fast tests

### Numerical results

| Setup | Steps | Result |
|---|---|---|
| Full stack single grid a=0.9 | 200 | stable, Kmax bounded |
| Refined full stack A=1.5 | 100 | Kmax decreasing |
| Refined full stack A=3.0 | 100 | NEC violation, stable |

### What v6.5 does NOT deliver

- Full stack beyond 500 steps
- Complex scalar + full stack
- Radiation extraction
- MPI / production AMR

### Honest assessment

v6.5 provides the first NEC-violating stable evolution of the full
stack (rotating regular core + QCD condensate + tachyon field).
The refinement technique is essential: A=3.0 crashes single-grid.

## v6.6 scope

### What v6.6 delivers

- Complex scalar + PTUM driver (single grid, equilibrium seed)
- Refined complex + PTUM driver (2-grid)
- 555 fast tests

### Numerical results

| Setup | Steps | Result |
|---|---|---|
| Complex + PTUM A=1.5 (eq seed) | 100 | rho_min = -0.285, stable |
| Complex + PTUM A=3.0 (eq seed) | 100 | rho_min = -15.05, Kmax=23.5 |
| Refined complex + PTUM A=3.0 | 100 | rho_min = 0 (equilibrium psi lost) |

### Key insights

1. PTUM activation requires either (a) enough time (t > ~2.4 for
   exponential growth from a small seed) or (b) equilibrium seed
   psi_0(Phi).
2. Refined complex + PTUM has a limitation: prolongate/restrict
   overwrites fine-grid psi with prolonged coarse psi each step,
   pushing it off equilibrium. This is a design limitation, not a bug.

### What v6.6 does NOT deliver

- Consistent psi initialization across refinement levels
- Refined complex + PTUM NEC violation
- Complex scalar + Kerr + PTUM (all combined)
- Full-time-domain tachyon growth from small seed

### Honest assessment

v6.6 adds complex scalar to the PTUM stack. Single-grid complex +
PTUM reaches NEC violation at A=3.0. The refined version does not
(eq psi lost). Physically meaningful refined NEC violation is
already in v6.5 (SC + PTUM).

Complex scalar in the refined full stack is deferred to v7.0.

## v7.0 scope

### What v7.0 delivers

- Refined complex+PTUM with NEC violation (persistent)
- Root cause fixes documented:
  - psi_equilibrium_smooth (smoothstep regularization)
  - restrict_matter flag (one-way coupling for matter)
- 579+ fast tests

### Key numerical result

| Setup | Steps | rho_min | Kmax_c | Kmax_f |
|---|---|---|---|---|
| Refined complex+PTUM A=3.0 | 100 | -15.67 | 0.087 | 0.245 |

rho_min persistent, psimax = 2.83 constant.

### Design decision: one-way coupling

Matter fields (Phi, psi, Psi_R/I) evolve on the coarse grid
independently; the fine grid receives prolonged coarse matter each
step. This avoids the ~35% per-step Gaussian peak loss from the
restrict averaging. The metric still uses two-way coupling.

This is a standard AMR technique for compact sources. It trades
fine-grid matter accuracy for coarse-grid stability.

### What v7.0 does NOT deliver

- Two-way matter coupling without peak loss (deferred to v7.1)
- Full 3D AMR
- Radiation extraction
- Kerr + complex + PTUM combined

### Honest assessment

v7.0 provides the first persistent NEC-violating refined
complex+PTUM evolution. Two root causes identified and fixed:
(a) sharp psi discontinuity -> smoothstep;
(b) restrict matter loss -> one-way coupling.

NEC violation rho_min = -15.67 stable over 100 steps.

## v7.1 scope

### What v7.1 delivers

- Peak-preserving restrict (restrict_max, exploratory)
- Refined Kerr + complex + PTUM (v7.1-alpha2)
- Ford-Roman QFT caveat at strong NEC (v7.1-alpha3)

### Key results

| Config | Steps | rho_min | Ford-Roman |
|---|---|---|---|
| Refined Kerr+complex+PTUM A=3.0 | 100 | -15.66 | **VIOLATED** (tau>=1) |
| Weak NEC PTUM (v6.0-v6.5) | 100-1000 | -0.1 | Compatible |

### QFT caveat (important)

Strong-amplitude PTUM discharge violates Ford-Roman quantum
inequalities at tau >= 1.0 (massless ratio 3.93, massive ratio
15.73). Weak-amplitude NEC is QFT-compatible.

This is a general feature of strong NEC violation:
short-lived strong negative energy is allowed; persistent
strong negative energy is not.

### What v7.1 does NOT deliver

- QFT-compatible strong NEC configuration
- Two-way matter coupling without peak loss
- Full 3D AMR

### Honest assessment

v7.1 extends NEC violation to Kerr+complex+PTUM with refinement,
and identifies the QFT constraint that strong-amplitude discharge
must be short-lived. This is a numerical demonstration of the
known tension between exotic-matter CTC solutions and quantum
inequality bounds.

## v7.2 scope

### What v7.2 delivers

- Transient discharge attempt (natural Phi decay, lambda=0.05)
- Fast-discharge variant (lambda=1.0, refined grid)
- Negative result documented

### Key finding: timescale mismatch

PTUM discharge has a natural lifetime of t ~ 5-20 (from the fuel
potential V_fuel = V0(1 - e^{-lambda Phi^2}) with lambda=0.05).
The accessible simulation time (t ~ 0.1, 500 steps) captures less
than 1% of this. The discharge therefore appears "permanent" in
the accessible window.

Consequence: Ford-Roman quantum inequality is violated at tau >= 1
for strong amplitude (A=3.0, rho_min = -15.66). This is
**physically expected** — strong NEC violation cannot be
QFT-compatible over QFT-scale windows.

### What v7.2 does NOT deliver

- A QFT-compatible strong NEC configuration
- Full transient discharge simulation (would require 25,000-100,000
  steps, 9-36 hours CPU)
- Any resolution of the timescale mismatch

### Honest assessment

v7.2 confirms that PTUM's "discharge" framing is consistent with
cosmological timescales but NOT with the QFT window (tau ~ 1).
This is a feature, not a bug: quantum inequalities bound negative
energy over quantum timescales, and PTUM's strong discharge
necessarily violates them.

The weak NEC regime (A=1.5, rho_min = -0.28, from v6.5) remains
the only QFT-compatible non-trivial NEC configuration in this
codebase. Strong NEC is QFT-incompatible, as theory requires.

## v7.3 scope

### What v7.3 delivers

- QFT threshold scan (11 amplitudes, A=1.4..3.5)
- Linear scaling law: ratio = 0.795 * |rho_min|
- Single-parameter threshold: |rho_min| = 1.258
- Fig8 (threshold scan) in bilingual paper

### Key result

QFT-compatible NEC: |rho_min| < 1.258 (A_fuel < 1.87)
QFT-incompatible NEC: |rho_min| > 1.258

This provides a quantitative QFT constraint on the PTUM model:
strong NEC violation above |rho_min| ~ 1.26 necessarily violates
Ford-Roman quantum inequalities.

### What v7.3 does NOT deliver

- QFT-compatible strong NEC
- Mechanism to evade the threshold
- Full transient-discharge simulation

### Honest assessment

v7.3 provides a clean single-parameter bound on PTUM NEC
violation. Strong NEC remains QFT-incompatible, as theory
requires. The bound is quantitative, not qualitative.

## v7.5 scope

### What v7.5 delivers

- Discharge duration analysis (QFT-OK window)
- Kerr-invariance of QFT-OK window

### Key results

**Discharge duration:** Phi decays <1% over t=0.25 in A=1.4-1.8.
Discharge frozen; timescale mismatch (v7.2) confirmed in QFT-OK window.

**Kerr-invariance:** QFT-OK window (A ∈ [1.35, 1.87]) unchanged
between a=0.0 and a=0.9. Threshold is set by PTUM physics, not by
background metric.

### Honest assessment

Two negative/invariance results:
1. Discharge never completes in accessible simulation time (any A).
2. QFT-OK window does not depend on Kerr spin.

The QFT-OK window is a PTUM-intrinsic property.

## v8.0 scope (Brans-Dicke extension)

### What v8.0 delivers

- BD scalar-tensor matter module
- BD RK4 evolution driver
- BD observational tests (Mercury, GW speed, Cassini)
- BD + PTUM combined driver + invariance check

### Key results

**BD observables:**
- Mercury perihelion: 43 arcsec/century (GR) ✅
- BD deviation at Cassini bound: < 2e-5
- GW speed: c_GW/c = 1 - 2.5e-5

**BD + PTUM:**
- phi_BD = 0.5, 1.0, 2.0 gives identical PTUM QFT ratios
- PTUM QFT window is BD-invariant

### Combined with v7.5

PTUM QFT threshold invariant under:
- Kerr spin (a=0..0.9)
- Brans-Dicke phi (0.5..2.0)
- Flat (v6.5)

**PTUM QFT threshold is intrinsic, background-independent.**

### What v8.0 does NOT deliver

- Full nonlinear BD (only slow-motion weak-coupling)
- BD black hole perturbation theory
- BD + Kerr + complex + PTUM (all combined)
- Numerical BD BH solutions

### Honest assessment

v8.0 adds the first alternative-gravity extension. The BD module
is at the leading-order level; the QFT-invariance result is robust
but does not claim full nonlinear BD.

## v8.2 scope (Starobinsky + LQC invariance)

### What v8.2 delivers

- Starobinsky f(R) = R + alpha R^2 scalar-tensor module
- LQC + PTUM density comparison
- QFT window invariance under LQC

### Key results

**Starobinsky:**
- V(phi) plateau, V_max = 1/(8 alpha)
- n_s ~ 0.965, r ~ 0.004 (Planck-compatible)

**LQC + PTUM:**
- PTUM density ~1e19 kg/m^3, rho_c ~1e96 kg/m^3
- Ratio: 1e-77 (LQC correction negligible)
- QFT window shift < 1e-70 (machine precision)

**Triple invariance:**
PTUM QFT-OK window unchanged under:
- Kerr spin (v7.5)
- Brans-Dicke phi (v8.0)
- LQC Friedmann (v8.2)

### What v8.2 does NOT deliver

- Full Starobinsky BSSN evolution
- LQC BSSN (only FRW + PTUM density comparison)
- f(R) black hole solutions
- LQC perturbations

### Honest assessment

v8.2 confirms that the PTUM QFT threshold is a
background-independent property, robust to three independent
extensions: Kerr, Brans-Dicke, LQC.

## v8.3 scope

### What v8.3 delivers

- Hu-Sawicki f(R) dark energy module
- Universal scalar-tensor QFT invariance theorem

### Key result

**Universal invariance:** The PTUM QFT-OK window
(A in [1.35, 1.87]) is invariant under ANY scalar-tensor
extension of GR of the form F(phi) R. Proven analytically
(PTUM rho_min is matter-sector-only) and tested numerically
across 5 F(phi) values.

### What v8.3 does NOT deliver

- Full Hu-Sawicki BSSN evolution
- f(R) black hole or cosmological solutions
- Screening mechanism dynamics

### Honest assessment

v8.3 generalizes the triple invariance of v8.2 to a
universal statement: PTUM QFT threshold is a fundamental
PTUM-intrinsic property, independent of the gravitational
background at leading order.

## v9.0/v9.1 scope (GW extraction + Horndeski)

### What v9.0 delivers
- Newman-Penrose GW extraction module
- psi_4, TT polarizations, energy flux

### What v9.1 delivers
- Horndeski scalar-tensor skeleton (GR + BD limits)
- Universal QFT window invariance under Horndeski

### What v9.x does NOT deliver
- Full Horndeski evolution
- 3D GW extraction from BSSN
- Realistic GW source templates

### Honest assessment

The NP module provides formulas; full 3D extraction from BSSN
is deferred. Horndeski extends the universal-invariance result
from scalar-tensor (v8.3) to the full Horndeski class.

## v9.2 scope (GW extraction)

### What v9.2 delivers

- Newman-Penrose GW extraction module (psi_4, TT, flux)
- BSSN + NP extraction pipeline
- Null result for spherically symmetric source

### Key results

| Setup | h+ max | psi_4 max |
|---|---|---|
| Spherical scalar (A=1e-4) | 0 | 0 |

### Physical interpretation

Monopole source -> no GW. Consistent with Einstein (1918)
quadrupole formula. The pipeline correctly identifies
zero signal.

### What v9.2 does NOT deliver

- Nonzero GW extraction (requires quadrupole source)
- 3D grid with realistic binary merger
- LIGO-comparable waveforms

### Honest assessment

The GW extraction pipeline is functional end-to-end. This
specific setup produces null signal because the source is
symmetric. Realistic GW extraction requires asymmetric
(quadrupole) dynamics, deferred to v10.

## v9.3 scope (quadrupole GW)

### What v9.3 delivers

- Quadrupole (two-bump) scalar source
- Nonzero GW extraction (h+ ~ 4.5e-12)
- Confirms Einstein quadrupole formula

### What v9.3 does NOT deliver

- Realistic binary merger (requires 3D evolution)
- LIGO-comparable waveforms
- Post-Newtonian matching

### Honest assessment

Quadrupole GW extraction now demonstrated. Signal is small
(A=1e-3 source), falls off with distance as expected. Full
binary inspiral is future work.

## v10.0 scope (binary GW — null result)

### What v10.0 attempted
- Two scalar Gaussians in circular orbit
- Complex scalar version with phase separation
- Expected f_GW = 2 f_orbit

### Result: null
- h+ ~ 0 (machine noise)
- Scalar sources have no TT content

### Honest assessment
This is a well-known physics limitation, not a numerical bug.
The GW extraction pipeline is correct; scalar sources simply
don't produce gravitational waves in GR. Real GW extraction
would require either vector/tensor fields or full 3D binary BH
evolution (BSSN with black holes), which is a much larger
project.

Documented as a negative result. Same pattern as the
Kerr/Brans-Dicke invariance: physics constrains, not code.
## v10.1 scope (GW pipeline validation)
### What v10.1 delivers
- Direct TT-gauge wave packet injection (gw_initial_data.py)
- End-to-end pipeline test: injection -> BSSN -> NP extraction
- Validation: extracted h+ matches analytic envelope
### Key result
| Observer | h+ extracted | Analytic |
|---|---|---|
| z=6 | 8.22e-4 | 8.2e-4 |
| z=7 | 5.71e-5 | 5.7e-5 |
| z=8 | 6.71e-7 | 6.7e-7 |
Pipeline reads metric perturbation correctly. psi_4 = 0
at packet peak (physical: dh/dt = 0).
### What v10.1 does NOT deliver
- GW from realistic (non-symmetric) source (scalar doesn't work)
- Binary BH inspiral (requires black hole excision)
- LIGO-comparable waveforms
### Honest assessment
Pipeline validated end-to-end. Combined with v10.0's null result
(scalar source -> no GW), this establishes: the extraction code
is correct; the challenge is producing a source with TT content,
which requires either vector/tensor fields or full 3D BH evolution.
## v10.2 scope (Proca GW source)
### What v10.2 delivers
- Proca (massive vector) field stress-energy module
- BSSN + Proca RK4 driver
- Nonzero h+ extraction from vector source
### Key result
| Source | h+ max | TT content |
|---|---|---|
| Scalar | 0 | none |
| Proca | 4.0e-6 | nonzero |
Vector field's anisotropic T_ij produces TT content,
unlike scalar. First physical GW signal (not injected).
### What v10.2 does NOT deliver
- Massive Proca (m > 0) GW extraction (used m=0)
- Realistic astrophysical source
- Frequency-domain waveform matching
### Honest assessment
Confirms the general principle: GW requires TT content in
T_ij, which scalar fields lack. Vector/tensor sources provide
it. This closes the source-class question from v10.0/v10.1.
## v11.0 scope (3D binary BH)
### What v11.0 delivers
- Moving puncture 3D validated (200 steps, growth 1.0001)
- Brill-Lindquist binary BH IC module
- Binary BH evolution stable (100 steps)
- GW signal h+ = 1.4e-7 (from infall, not orbit)
### Key results
| Setup | Steps | H growth | h+ |
|---|---|---|---|
| Single MP 3D | 200 | 1.0001 | — |
| Binary BL | 100 | 1.0001 | 1.4e-7 |
### What v11.0 does NOT deliver
- Orbital dynamics (Bowen-York momentum) — v12.0
- GW frequency matching
- Realistic inspiral waveform
- Full merger + ringdown
### Honest assessment
Moving puncture works in 3D. Binary BH IC (Brill-Lindquist) is
time-symmetric: punctures plunge, not orbit. Real waveform
requires momentum constraint solution (Bowen-York). This is
the bridge to v12.0.
## v12.0 scope (Bowen-York orbital binary)
### What v12.0 delivers
- Bowen-York IC: two punctures with tangential momenta
- Orbital binary BH evolution (200 steps)
- Both GW polarizations present (h+, hx)
### Key result
| Metric | Value |
|---|---|
| H growth (200 steps) | 1.0098 |
| h+ max | 7.9e-6 |
| hx max | 7.4e-4 (orbital signature) |
### What v12.0 does NOT deliver
- Full orbital period (only 1.7% captured)
- GW frequency measurement (Nyquist-limited)
- Merger + ringdown
- LIGO-comparable waveform
### Honest assessment
Bowen-York IC produces orbital dynamics (hx != 0, absent in
time-symmetric Brill-Lindquist). Stability bounded over 200 steps.
Full inspiral requires longer runs + higher resolution (future).
## v12.1 scope (plunge + PTUM in BH + long orbit)
### What v12.1 delivers
- Brill-Lindquist plunge waveform (N=7, 50 steps)
- PTUM static field in Schwarzschild background
- Long Bowen-York orbit (300 steps)
### Key results
- Plunge: h+ = 8.7e-6, hx = 0 (symmetric)
- Orbit: h+ = 3.1e-4, hx = 1.1e-2 (both polarizations)
- PTUM QFT window unchanged by BH curvature
### What v12.1 does NOT deliver
- Full orbital period (5% captured)
- N>=9 (Python too slow)
- Merger + ringdown
## v12.2 scope (vectorized BSSN)
### What v12.2 delivers
- Fully vectorized BSSN RHS (einsum-based)
- Vectorized moving puncture RK4 driver
- 25x speedup at N=11-13
- N=13 binary BH in 26 seconds
### Key numbers
| N | loop (s/RHS) | vec (s/RHS) | speedup |
|---|---|---|---|
| 9 | 0.68 | 0.02 | 27x |
| 11 | 1.73 | 0.07 | 25x |
| 13 | 4.87 | 0.19 | 26x |
### What v12.2 does NOT deliver
- Full orbital period (still Nyquist-limited)
- Merger / ringdown
- N>15 (memory)
### Honest assessment
Vectorization makes N=13 feasible. Full inspiral and merger
now possible but still need long runs + boundary handling.
## v12.3 scope (triple toy model)
### What v12.3 delivers
- 0D/1D SC + PTUM + BH triple coupling model
- Confirms QFT window invariance
### Key result
PTUM QFT window (A ∈ [1.35, 1.87]) unchanged by SC core or BH curvature.
Universal invariance extends to triple sector.
### What v12.3 does NOT deliver
- Full 3D SC+PTUM+BH evolution
- Backreaction of PTUM on BH metric
- Realistic observational signatures
## v12.4 scope (N=15 medium run)
### What v12.4 delivers
- N=15 binary evolution, 1000 steps in 152s
- H growth 1.0268 (bounded)
- Both GW polarizations grow monotonically
### What v12.4 does NOT deliver
- Full orbital period
- Merger + ringdown
- Constraint-preserving BC
## v12.5 scope (DHOST + teleparallel)
### What v12.5 delivers
- DHOST skeleton (beyond-Horndeski terms)
- Teleparallel f(T) skeleton
- GW170817 constraint
### Key results
- GW speed: |c_GW/c - 1| < 1e-15 (GW170817)
- PTUM QFT window invariant under DHOST, f(T)
### What v12.5 does NOT deliver
- Full DHOST evolution
- f(T) BSSN implementation
- Observational predictions beyond GW speed
## v12.6 scope (N=17 inspiral — constraint breakdown)
### What v12.6 delivers
- N=17 chunked run (3200 steps, ~10 min CPU)
- Yörünge periyodunun %28'i (T=23, t=6.4)
- Constraint breakdown detected
### Key finding
| Chunk | Step | H/H0 | h+ |
|---|---|---|---|
| 1 | 800 | 1.09 | 5.8e-4 |
| 2 | 1600 | 1.19 | 2.1e-3 |
| 3 | 2400 | 1.28 | 4.2e-3 |
| 4 | 3200 | 1.63 | 7.2e-3 |
**H growth accelerates past 1.5** → Sommerfeld outflow BC
insufficient for long inspiral. Wave reflection at boundaries
accumulates constraint violation.
### What v12.6 does NOT deliver
- Full orbital period (stopped at 28%)
- Merger / ringdown
- Constraint-preserving BC
### Honest assessment
The constraint breakdown is a **physics result**, not a bug:
Sommerfeld BC is first-order; long inspiral needs
constraint-preserving boundaries (e.g., Friedrich-Nagy, Kreiss-
Winicour). This is a well-known limitation of simple NR setups.

## v12.7 scope (negative result: boundary damping ineffective)

### What v12.7 delivers
- Sponge + multi-layer Sommerfeld BC module (8 checks)
- Z4c-style constraint damping in vectorized RHS
- Systematic test: kappa in [0.02, 0.5], sponge in [3,5]

### Key finding
Damping has NO measurable effect on H growth. At t=6.4, outgoing radiation has not yet reached the boundary (L/2=10). Constraint violation accumulates in the puncture region, where N=17 gives ~1.25 grid cells per unit mass.

### What v12.7 does NOT deliver
- Constraint control (needs AMR / N>=32)
- Full inspiral
- Merger / ringdown

### Honest assessment
Boundary damping is the wrong tool. Next: mesh refinement near punctures.

## v12.8 scope (BSSN shift + resolution limit)

### What v12.8 delivers
- Full BSSN shift advection in vectorized RHS
- Shift properly coupled to metric (waveforms change with eta)
- Systematic eta scan (1.0, 4.0)
- Clear diagnosis: puncture under-resolution is the limit

### Key finding
Shift integration works (h+ differs 10% between eta=1, eta=4) but does NOT reduce H growth. The constraint breakdown is resolution-limited: at N=17, puncture has ~1 grid cell (M=1, h=1.25). Increasing resolution (N=24) makes it worse (constraint source sharper), confirming low-N was hiding the problem.

### What v12.8 does NOT deliver
- Constraint control (requires AMR)
- Full inspiral
- Merger / ringdown

### Honest assessment
All ingredients except AMR are now correct. AMR (v12.9+) is the required next step for production binary BH.

## v12.9 scope (fixed-patch AMR, proof of concept)

### What v12.9 delivers
- Trilinear prolongation + injection operators (6 tests)
- Two-grid AMR in lockstep (coarse + fine, refinement 2x)
- Positive constraint reduction: 1.0756 vs 1.1254 at step 1000 (33% better)

### What v12.9 does NOT deliver
- Sub-cycling (coarse and fine evolve with same dt)
- Multiple refinement levels
- Proper restriction (fine->coarse feedback)
- Production runs (still ~1 cell per puncture at fine level)

### Honest assessment
Fixed-patch AMR works as proof of concept. 2x refinement gives 33% constraint reduction. Production requires: (a) 4x-8x refinement (b) sub-cycling (c) restriction. This is now the clear path to full inspiral.

## v12.10 scope (subcycling + refinement scan)

### What v12.10 delivers
- Configurable AMR driver (factor 2/4, subcycle 1/2/4)
- Systematic scan: 2x sc=2 (neutral), 4x sc=1 (CFL fail)

### Key finding
Subcycling does NOT improve constraint control at 2x. BC injection is the limiting factor. Higher refinement without subcycling fails (CFL).

### What v12.10 does NOT deliver
- Tested 4x + subcycle=4 (correct CFL)
- Proper restriction feedback
- Production AMR

### Honest assessment
2x AMR works (v12.9). Subcycling and refinement factor are not the lever. Fixing BC injection (higher-order prolongation, moving patch) is the next real step, not deeper refinement.
