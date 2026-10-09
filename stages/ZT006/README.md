# ZT-006 — Stress-Energy and Linearized-GR Source Chain

Stage goal: build the electromagnetic stress-energy tensor from validated
fields, condition it numerically, and couple it to a linearized-GR Green
kernel with explicit convergence diagnostics.

## Honest status of record

* **ZT-006.2.19 — FAIL (kept on record).** The final-grid Green-kernel
  result was NOT converged for observers inside the source volume
  (0/20 grid checks; last-grid changes of ~4.8% G1, ~37% G2, ~80% G3,
  ~25% G4). Nothing about this stage may be labeled "successful".
* **ZT-006.2.20 — FAIL-equivalent (kept on record).** Moving all
  observers outside the source volume removed the 1/r singularity but did
  NOT fix convergence. The defect therefore is not merely the interior
  singularity; suspicion moved to the source distribution and/or spatial
  discretization / field differentiation. See
  `results/archived/zt006_2_20_results.json`.

## Consequences (implemented in ZT-006.3)

The frozen/instantaneous Green kernel is not acceptable as the source-to-
metric transition. ZT-006.3 adds:

* a retarded linearized-Einstein source solver
  (`src/zt005/retarded_gr.py`),
* a stress-energy conservation diagnostic ∂_μ T^{μν}
  (`src/zt005/stress_energy_conservation.py`),
* a Lorenz-gauge consistency diagnostic ∂_μ h̄^{μν}
  (`src/zt005/gauge_checks.py`),
* known-metric / geodesic / CTC benchmark layers
  (`src/zt005/benchmark_metrics.py`, `geodesic_solver.py`,
  `causal_structure.py`, `ctc_search.py`).

FAIL results are data. They are archived under `results/archived/` and
never deleted.

## ZT-006.2.24 — RESOLVED

External-observer Green-kernel convergence now passes.

- Root cause: filament source at 5 mm radius left T(x') too sharp for
  the quadrature; consecutive-grid rel error stayed O(0.1-0.5).
- Fix: effective smoothing radius `a_wire_eff = 0.08 m` in the
  integration kernel (numerical, not physical wire radius).
- Result: `rel < 0.05` on all consecutive grid pairs, all four
  source configs (G1-G4). Final step rel ~ 5.8e-06 (full mode).
- Evidence: `results/validated/zt006_2_24_results.json`
- Script: `src/zt005/run_zt006_2_24.py`
