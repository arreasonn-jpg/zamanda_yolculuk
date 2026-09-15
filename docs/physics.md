# Physics model and conventions

## What this project is (and is not)

This repository is a research codebase for **numerically testing a
spacetime-manipulation hypothesis**: can a driven electromagnetic device
produce a measurable, causally interesting metric perturbation?

It is **not** a time machine, a CTC generator, or a device design. The
wording "time machine" is deliberately avoided in technical layers; the
correct description of the current work is *spacetime-manipulation
hypothesis testing*.

## Layer structure

1. **Electromagnetic layer (ZT-005)** — magnetoquasistatic (MQS) model:
   vector potential of closed filamentary currents, harmonic drive,
   internal consistency checks (div B ≈ 0, curl E ≈ −dB/dt).
2. **Stress-energy layer (ZT-006)** — T^{μν} from the EM fields only
   (see assumptions.md: this is an *incomplete source model*).
3. **Linearized-GR layer (ZT-006.3)** — retarded solution of the
   linearized Einstein equation on a fixed Minkowski background:

       □ h̄_μν = −16πG/c⁴ T_μν ,   ∂_μ h̄^{μν} = 0 (Lorenz gauge)

       h̄_μν(t, x) = (4G/c⁴) ∫ T_μν(t − |x−x′|/c, x′) / |x−x′| d³x′

   The integrand is **T_μν, not its trace reverse**; trace reversal lives
   on the field side (h̄_μν = h_μν − ½ η_μν h).
4. **Benchmark/geometry layer (ZT-006.3)** — known exact metrics,
   geodesic integration, causal classification, closed-timelike-curve
   known-answer tests. This layer validates the machinery before it is
   ever pointed at device-generated fields.

## Conventions

| Item | Legacy EM layer (ZT-005/6) | New GR layer (ZT-006.3+) |
|---|---|---|
| Signature | (+, −, −, −) | (−, +, +, +) |
| Timelike interval | ds² > 0 | **ds² < 0** |
| Coordinates | SI, x⁰ = ct implicit | x^μ = (ct, x, y, z), SI |
| Units | SI | SI for the retarded solver; geometric units c = G = 1 for benchmark metrics |

Trace-reversed tensor *components* are identical in both signature
conventions, so stress-energy arrays flow between layers unchanged.

## Key equations of record

* Retarded linearized solution (BLOCKER item of the ZT-006.2 postmortem):

      h̄_μν(t,x) = (4G/c⁴) ∫ T_μν(t−|x−x′|/c, x′)/|x−x′| d³x′

* Conservation requirement for a trustworthy Einstein source:

      ∂_μ T^{μν} = 0     (numerically diagnosed in
                           stress_energy_conservation.py)

* Gauge condition:

      ∂_μ h̄^{μν} = 0     (numerically diagnosed in gauge_checks.py)

* Geodesic equation:

      d²x^μ/dτ² + Γ^μ_{αβ} (dx^α/dτ)(dx^β/dτ) = 0

* CTC known-answer thresholds:
  * Gödel universe: φ-circles timelike for r > asinh(1) = ln(1+√2).
  * Idealized Tipler cylinder: φ-circles timelike for a·r > 1.

## Benchmark metrics

| Metric | Coordinates | CTC status | Role |
|---|---|---|---|
| Minkowski | Cartesian | none | flat reference |
| Schwarzschild | (t, r, θ, φ), r > 2M | none exterior | precession benchmark |
| Kerr (Boyer–Lindquist) | (t, r, θ, φ), r > r₊ | none exterior | rotation benchmark |
| Gödel (cylindrical) | (t, r, φ, z) | yes, r > asinh(1) | positive CTC answer |
| Tipler cylinder (idealized) | (t, r, φ, z) | yes, ar > 1 | positive CTC answer |
| Morris–Thorne (static) | (t, l, θ, φ) | none in static form | topology benchmark |

The Tipler entry is a pedagogical idealization (infinite rotating
cylinder) used strictly as a known-answer test for the detector.

## What the h ~ 10⁻⁴⁶ number means

Published perturbation magnitudes of order h ~ 10⁻⁴⁶ under the current
normalization are a physics *result*, not a bug: they quantify how far
the device is from producing any macroscopic metric effect. A final
project conclusion of the form *"mathematically consistent, but the
required energy scale is practically inaccessible"* is a valuable,
publishable outcome — not a failure.
