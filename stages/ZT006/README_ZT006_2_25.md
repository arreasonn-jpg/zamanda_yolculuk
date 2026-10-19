# ZT-006.2.25 — Retarded-Solver Device-Source Convergence

Closes the open item left in `docs/validation.md`:
grid-refinement study of the **retarded** integral with the actual
device source distribution (the follow-up to ZT-006.2.19/.20,
previously resolved only for the instantaneous kernel in
ZT-006.2.24).

## Run

```bash
python -m zt005.run_zt006_2_25 --mode quick
python -m zt005.run_zt006_2_25 --mode full
```

## Setup

* Solver: `zt005.retarded_gr.retarded_hbar` (trace-reversed,
  (-+++) layer, singular-cell correction).
* Sources: device G1–G4 (`geometry.build_sources`), external
  observers (same set as ZT-006.2.24), quadrature from
  `focused_source_nodes`, effective smoothing `a_wire_eff = 0.08 m`
  (same as 2.24).
* Time dependence (MQS 1 MHz drive):
  `A(x,t) = A1 cos(ωt_ret)`, `B(x,t) = B1 cos(ωt_ret)`,
  `E(x,t) = ω A1 sin(ωt_ret)` with per-point retarded time
  `t_ret = t − R/c`. Across the source ωL/c ≈ 5e−8 rad — formally
  kept, far below quadrature error.
* Observation times: half cycle at 1 MHz (0, T/16, T/8, T/4,
  3T/8, T/2).

## Result (quick mode, 2026-10-19)

* All four source configs (G1–G4): rel < 0.05 on consecutive
  grid pairs; worst rel ≈ 1.47e−2, best ≈ 3.47e−3.
  **CONVERGENCE: PASS.**
* Retarded residuals are numerically indistinguishable from the
  instantaneous ZT-006.2.24 values (as expected for MQS:
  λ ≈ 300 m ≫ L ≈ 2.5 m).

## Honest caveats

* This validates the *spatial* convergence of the retarded
  integral with a regularized device source. It does NOT
  constitute a validated device *metric* claim: T_μν is still
  EM-only (incomplete source model) and gravity remains
  linearized; expectation stays h ~ 1e−46 scale.
* The a_wire_eff = 0.08 m smoothing is a numerical regularizer,
  not a physical wire radius (unchanged from 2.24).

Artifact: `results/validated/zt006_2_25_results.json`.
