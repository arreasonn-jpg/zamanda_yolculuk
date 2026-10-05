# Phase plan (locked)

Progression is gated: a phase starts only when the previous phase's
validation exists on record.

| Phase | Content | Status |
|---|---|---|
| 0 | Project infrastructure (repo layout, CI, docs, results policy) | ✅ this reorganization |
| 1 | Maxwell validation (MQS consistency) | ✅ ZT-005.x (full-wave deferred to ZT-007) |
| 2 | T_μν validation | 🔶 partial: EM-only T built; conservation diagnostic now exists but device-source convergence FAILED (ZT-006.2.19/.20) |
| 3 | Retarded linearized Einstein solver | 🔶 solver implemented + manufactured-source validated; device-source convergence study pending |
| 4 | Metric validation | ✅ benchmark metrics validated (known metrics) |
| 5 | Geodesic engine | ✅ RK4 integrator validated (norm, precession) |
| 6 | Known-CTC benchmark reproduction | ✅ Gödel + Tipler known answers reproduced; negative benchmarks clean |
| 7 | Device-generated spacetime search | ❌ not started — gated on Phases 2/3 convergence |
| 8 | Energy scaling (P_in/loss/stored/radiated/thermal, η, E_{1s}) | ❌ not started |
| 9 | Earth + target-coordinate mapping (tensor level) | ❌ kinematic placeholder only |
| 10 | Engineering feasibility (coils, conductors, thermal, HV, power, backpack constraint optimizer) | ❌ not started |
| 11 | Experimental apparatus (instrument-only; no human subjects) | ❌ not started |

## Immediate next moves (in order)

1. **Conservation of the real device source.** Evaluate ∂_μ T^{μν} for
   the MQS device fields (including the wire currents explicitly) and
   make the residual a tracked regression.
2. **Retarded-solver convergence for the device source.** The failure
   mode of ZT-006.2 must not repeat: grid-refinement study with
   analytic/semi-analytic sub-checks before any metric claim.
3. **Gauge audit of the device pipeline.** ∂_μ h̄^{μν} for the device
   retarded solution.
4. **Energy bookkeeping skeleton** (Phase 8 constants and budgets), so
   that the h ~ 10⁻⁴⁶ scale can be translated into an explicit energy
   requirement.

## Explicitly deferred

* C/C++ kernels and GPU acceleration — not before the Python baseline is
  right (a fast wrong answer is worse than a slow right one).
* Docker / packaging polish — after Phase 3 closure.
* Full-wave Maxwell (ZT-007) — when frequency scans demand it.

## Human-experiment policy

Any human experiment is out of scope for all listed phases and remains
so until an independent ethics/safety review concludes otherwise.
