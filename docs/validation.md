# Validation status

Rule: FAIL results are preserved and cited. Nothing is labeled
"successful" unless its acceptance test passes.

## Current state (ZT-006.3)

### PASS — EM layer (ZT-005.x)

MQS consistency suite: div B ≈ 0, curl E ≈ −dB/dt, curl A vs Biot–Savart,
filament/grid convergence for the EM fields. See
`results/archived/zt005_*.json` and `stages/ZT005/`.

### FAIL — frozen Green kernel convergence (ZT-006.2.19 / .20)

* ZT-006.2.19: **0/20** grid-convergence checks passed with observers
  inside the source volume (last-grid changes ≈ 4.8% / 37% / 80% / 25%
  across G1–G4).
* ZT-006.2.20: external observers removed the 1/r singularity but
  convergence was still **not** achieved → the defect is not the interior
  singularity; suspicion shifted to the source distribution and/or
  spatial discretization / field differentiation.
* Consequence: the frozen kernel is disqualified as the source-to-metric
  transition; the retarded solver (below) replaces it.

### PASS — GR validation layer (ZT-006.3), 16/16 checks

Regenerate: `python -m zt005.run_zt006_3_gr_layer --mode quick|full`
(artifact: `results/validated/zt006_3_gr_layer_results.json`).

| Check group | Content | Status |
|---|---|---|
| Metric sanity | symmetry + Lorentzian signature of all 6 benchmarks | PASS |
| Geodesic norm conservation | g(u,u) drift < 1e−5 (numeric Γ), < 1e−9 (analytic) | PASS |
| Schwarzschild perihelion | measured advance vs 6πM/(a(1−e²)), < 5% | PASS |
| CTC known-answer (Gödel) | detection flips exactly at r = asinh(1) | PASS |
| CTC known-answer (Tipler) | detection flips exactly at ar = 1 | PASS |
| CTC negative benchmarks | no CTC found in Minkowski/Schwarzschild/MT | PASS |
| Conservation diagnostic | conserved manufactured source → small residual; broken source → large | PASS |
| Lorenz-gauge diagnostic | retarded solution of conserved source → ∂_μ h̄^{μν} small | PASS |

### PASS — retarded solver unit properties

* Causality: pulse observed only after t₀ + R/c, peak at the retarded
  time.
* Static limit: time-independent source reproduces the instantaneous
  integral to 1e−9.

## Not yet validated (open)

* Convergence of the retarded integral for the *device* source
  distribution (grid refinement study pending — this is where ZT-006.2
  died once; it must be redone against the retarded solver).
* Any nonlinear / self-consistent metric statement — none exists yet.
* Any CTC statement about device fields — the detector is validated only
  on known metrics so far.
* Energy model (P_in, P_loss, P_stored, P_radiated, P_thermal, η,
  E_{1s}) — not implemented.
