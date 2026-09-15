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
