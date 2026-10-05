# ZAMANDA YOLCULUK — spacetime-manipulation hypothesis testing

A research codebase that tests, numerically and honestly, whether a
driven electromagnetic device can produce a measurable metric
perturbation — and how large the gap to any causally interesting effect
really is.

> **Status in one line:** the numerical-validation core is strong; the
> physical source→metric chain is linearized, EM-only, and not yet
> converged for the device source. No time machine has been built, and
> the honest expectation is a quantitative bound, not a device.

The term "time machine" is avoided in technical layers by design: the
correct description is **spacetime-manipulation hypothesis**. All CTC
machinery in this repo is a known-answer benchmark, not a device claim.

## Repository layout

```
├── src/zt005/          # all code (EM layer, T_mu_nu, retarded GR layer,
│                       # benchmark metrics, geodesics, CTC search)
├── tests/              # pytest suite
├── stages/             # per-stage records (ZT001…ZT007), incl. FAILs
├── results/
│   ├── validated/      # acceptance-tested outputs (regenerable)
│   ├── exploratory/    # in-progress, not citable
│   └── archived/       # historical outputs — FAILs kept forever
├── docs/               # physics, assumptions, numerics, validation,
│                       # limitations, roadmap
├── configs/            # device baseline parameters
├── scripts/            # operational notes
└── .github/workflows/  # CI: pytest on push/PR
```

## Physics layers

| Layer | Module(s) | Status |
|---|---|---|
| Magnetoquasistatic EM (1 MHz, L ≈ 2.5 m, λ ≈ 300 m) | `maxwell`, `model`, `finite_wire_*` | validated (ZT-005.x) |
| EM stress-energy T_μν (**incomplete source model**) | `validated_stress_energy` | built; device-source convergence FAILED once (see below) |
| Retarded linearized Einstein solver h̄_μν(t, x) | `retarded_gr` | implemented; manufactured-source validated |
| Conservation ∂_μT^{μν} / Lorenz gauge ∂_μh̄^{μν} | `stress_energy_conservation`, `gauge_checks` | implemented with manufactured-source tests |
| Known metrics (Minkowski, Schwarzschild, Kerr, Gödel, Tipler, Morris–Thorne) | `benchmark_metrics` | validated |
| Geodesic integrator | `geodesic_solver` | validated (norm drift, perihelion advance) |
| Causal classifier + CTC known-answer search | `causal_structure`, `ctc_search` | validated on known metrics only |

## Honesty record (kept, not hidden)

* **ZT-006.2.19 — FAIL.** Green-kernel grid convergence: 0/20 checks
  with interior observers (last-grid changes ~4.8% / 37% / 80% / 25%).
* **ZT-006.2.20 — not converged.** External observers removed the 1/r
  singularity but convergence still failed → defect is in the source
  distribution and/or discretization, not only the singularity.
  Record: `results/archived/zt006_2_20_results.json`,
  `stages/ZT006/README.md`.
* Current perturbation scale h ~ 10⁻⁴⁶: a physics result quantifying the
  gap, not a bug. A final verdict of "consistent model, inaccessible
  energy scale" would be a valuable outcome.

See `docs/validation.md` for the full PASS/FAIL register and
`docs/assumptions.md` for every load-bearing assumption (EM-only source,
MQS regime, linearized gravity, resonance-as-hypothesis, kinematic Earth
model).

## Quick start

```bash
pip install -e .
pytest                                          # full suite
python -m zt005.run_zt006_3_gr_layer --mode quick   # GR validation layer
```

Key artifacts:

* `results/validated/zt006_3_gr_layer_results.json` — 16-check GR layer
  report (metric sanity, geodesic norm conservation, Schwarzschild
  perihelion vs 1PN, Gödel/Tipler CTC known answers, conservation and
  gauge residuals).

## Conventions

* New GR layer: signature (−,+,+,+), timelike = ds² < 0, x^μ = (ct, x, y, z).
* Legacy EM layer: (+,−,−,−). Trace-reversed components are identical in
  both, so arrays flow between layers unchanged. Details:
  `docs/physics.md`.

## Phase plan

Locked progression: infrastructure → Maxwell → T_μν → **retarded
linearized Einstein** → metric validation → geodesics → known-CTC
reproduction → device spacetime search → energy scaling → Earth/target
mapping → engineering feasibility → instruments. Human experiments are
out of scope at every phase. Full table: `docs/roadmap.md`.
