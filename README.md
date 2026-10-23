# ZAMANDA YOLCULUK — spacetime-manipulation hypothesis testing

A research codebase that tests, numerically and honestly, whether a
driven electromagnetic device can produce a measurable metric
perturbation — and how large the gap to any causally interesting effect
really is.

> **Status in one line:** the numerical-validation core is strong; the
> physical source→metric chain is linearized, EM-only, and not yet
> converged for the device source. No time machine has been built, and
> the honest expectation is a quantitative bound, not a device.
>
> **Current frontline (v14.0.0):** (a) device-source convergence is
> CLOSED — ZT-006.2.24 (instantaneous kernel) and ZT-006.2.25 (retarded
> solver, device sources) both PASS; the remaining gap to a device
> metric is physical (h ~ 10⁻⁴⁶, EM-only linearized chain), not
> numerical; (b) fixed-patch AMR for BSSN binaries: cubic prolongation
> (v12.11), finer fine resolution (v12.13), Gamma rederivation (v12.15)
> are all neutral/negative; COM-following moving patch quantified as a
> no-op (v13.0.0: COM drift < 1 coarse cell in 1000 steps; H_f/H0 =
> 1.0756 unchanged); 4x + subcycle=4 closed as NEGATIVE (v13.1.0:
> H@500 = 1.1418 vs baseline 1.0166 — refinement amplifies the
> prolonged-boundary mismatch); first POSITIVE lever is constraint
> projection (v13.2.0): Jacobi conformal solve of the injected phi
> band reduces H_f/H0 to 0.9319 @ 1000 (baseline 1.0756, -13.4%),
> below the run's own step-0 reference. Covariant-operator upgrade
> (v13.3.0) is NEUTRAL (0.9382 vs 0.9319): the ~13% gain is robust
> to the operator; iso-reference splits the sub-1.0 effect into
> ~10% definitional + rest genuine. Momentum branch (v13.4.0) also
> NEUTRAL: the band A-divergence is cleanable (per-step ratio 0.43)
> but carries no H signal — the phi-band projection is the whole
> effect; AMR lever list exhausted.
> (b2) full-inspiral stage-0 (v14.0.0): the first instrumented long
> run (no projection) blows up at t = 5.6 — the chronic AMR injection
> defect goes NONLINEAR in the injection BAND (fine band max |H|:
> 0.44 t=1 → 2.32 t=5, e-fold ~0.4 M; coarse grid stays healthy).
> The orbit-gauge signal itself is real (phase +0.44 rad, omega
> +0.13..0.15 rad/M ≈ Newtonian 0.125) but the band is the blocker:
> the measured critical path is projection-on long runs (v14.1).
> (c) code layout: BSSN-amr machinery in `zt005.amr/`, matter models
> in `zt005.matter/`.

The term "time machine" is avoided in technical layers by design: the
correct description is **spacetime-manipulation hypothesis**. All CTC
machinery in this repo is a known-answer benchmark, not a device claim.

## Version status

| Version | Highlights |
|---|---|
| v1.0.0 | Linearized GR, exotic matter models, energy gap report (43 orders) |
| v2.0.0 | Full BSSN solver (Gamma + alpha + beta + gauge + Z4c-style damping + physical Ricci); Schwarzschild/Kerr benchmarks stable |
| v3–v5   | Moving puncture, full Kerr (KS, a≤0.99), SC regular core, scalar matter, PTUM tachyon coupling (NEC violation confirmed) |
| v6–v8   | 2-level refinement (AMR-lite), Ford-Roman QFT bounds, QFT-OK window (A ∈ [1.35, 1.87]), Kerr/BD/LQC invariances |
| v9–v10  | Newman-Penrose GW extraction; scalar-source null result; end-to-end GW pipeline validation |
| v11–v12 | 3D binary BH: Brill-Lindquist → Bowen-York orbital; vectorized RHS (25x); N=17 constraint breakdown diagnosed; AMR v12.9–v12.13 (proof of concept, negative results documented) |
| v12.14  | ZT-006.2.25: retarded-solver device-source convergence PASS (source chain numerically closed) |
| v12.15  | Constraint-consistent Gamma injection: neutral (H unchanged); 737 fast + 21 slow tests; frontline: constraint projection / moving-patch AMR |
| v13.0.0 | COM-tracked moving patch quantified as no-op (COM drift < 1 coarse cell in 1000 steps; H_f/H0 = 1.0756 unchanged); subpackage restructure (`zt005.amr/`, `zt005.matter/`) |
| v13.1.0 | 4×+subcycle=4 (CFL-correct) closed as negative: H@500 = 1.1418 vs baseline 1.0166 — resolution amplifies the prolonged-boundary mismatch |
| v13.2.0 | First positive lever: band-local conformal projection (Jacobi solve of the injected phi band, guard-masked) — H_f/H0 @ 1000 = 0.9319 vs baseline 1.0756 (−13.4%), below the run's own step-0 reference |
| v13.3.0 | Exact γ̃-covariant operator + iso-reference metric: NEUTRAL (0.9382 vs 0.9319 flat) — the projection gain is operator-robust; ~10% of the sub-1.0 effect is definitional |
| v13.4.0 | Momentum Vector-Poisson band correction: NEUTRAL on H (0.9371 vs 0.9382) — band Ã-divergence cleanable (ratio 0.43/step) but carries no H signal; AMR lever longlist exhausted |
| v14.0.0 | Full-inspiral stage-0: unprojected long run blows up at t=5.6 IN THE BAND (fine band max \|H\| 0.44→2.32, e-fold ~0.4 M; coarse healthy); orbit-gauge signal real (ω +0.13..0.15 ≈ Newtonian); critical path = projection-on long runs |

See `docs/limitations.md` for the full scope/NOT-delivered register per
version, and `docs/validation.md` for the PASS/FAIL register.

## Repository layout

```
├── src/zt005/          # code root: EM layer, T_mu_nu, retarded GR layer,
│   │                   # benchmark metrics, geodesics, CTC search,
│   │                   # BSSN solver + run_* stage entry points
│   ├── zt005/amr/      # AMR-lite: interpolation/prolongation, BH tracking,
│   │                   # refinement, consistent injection, moving window
│   └── zt005/matter/   # matter models: scalar/complex/ptum/proca/
│                       # Brans-Dicke, condensate, exotic matter
├── tests/              # pytest suite
├── stages/             # per-stage records (ZT001…ZT010), incl. FAILs
├── paper/              # LaTeX sources + figures (arxiv_v3 active;
│                       # archive/ holds historical submission zips)
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
* **ZT-006.2.24/.25 — RESOLVED.** Effective smoothing `a_wire_eff`
  = 0.08 m recovers external-observer convergence (2.24), and the
  retarded solver passes the same device-source grid study (2.25,
  worst rel ≈ 1.5e−2 quick / 5.8e−6 full). Numerical chain closed;
  the gap to any device metric is physical, not numerical.
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
