# Assumptions register

Every assumption below is load-bearing. When one is removed, the
corresponding limitation in limitations.md is removed too.

## A1 — Incomplete source model (EM only)  ⚠️ HARD LIMIT

The current stress-energy tensor is

    T_μν = T_μν^EM

only. A physical device requires at least

    T_μν^total = T_μν^EM + T_μν^matter + T_μν^mechanical
                 + T_μν^thermal + T_μν^control + …

coils, conductors, structural stresses, thermal energy, the power source,
and pressure/stress contributions all gravitate. Until this is addressed
the source model is **incomplete by construction** and every metric
statement inherits that incompleteness.

## A2 — Magnetoquasistatic EM (not full-wave)

ZT-005 is a controlled MQS model. At f = 1 MHz and L ≈ 2.5 m,
λ ≈ 300 m, so L/λ ≈ 0.008 and MQS is defensible. Frequency scans that
raise L/λ require the full-wave stage (ZT-007).

## A3 — Linearized gravity on fixed Minkowski background

The retarded solver computes h̄_μν on flat space. There is no nonlinear
Einstein solve, no self-consistent g_μν(t, x) field yet, and no metric
back-reaction. All statements are first order in G.

## A4 — Resonance enhancement is a HYPOTHESIS, not a result

"Resonance strengthens the spacetime coupling" is an open hypothesis.
There is currently no model of Q, ω₀, L, C, R, or of the electromagnetic
modal structure. Code and documents must label resonance effects as
hypothetical until a modal model exists.

## A5 — Earth model is a kinematic placeholder

`target_time_mapping.py` provides a reference-frame-only model: circular
1 AU orbit, uniform sidereal rotation, spherical Earth, no tilt. It
contains **no tensor-gravity contribution** of the Earth, Sun, or
environment. The eventual decomposition

    T_μν^total = T_μν^device + T_μν^Earth + T_μν^solar + T_μν^environment

and the matching metric matching are Phase 9 items, not implemented.

## A6 — Device geometry is idealized

Filamentary closed loops on an idealized coil/backpack geometry. Finite
conductor cross-sections, real winding packs, and insulation are
engineering items (Phase 10), absent from the physics layer.

## A7 — Benchmark metrics are exact inputs, not device claims

Gödel / Tipler / wormhole metrics in the code exist only as known-answer
tests. Finding a CTC there validates the detector; it implies nothing
about the device.

## A8 — No human experiment

No stage of this project plans or permits a human subject. Any future
experimental apparatus phase is instrument-only and requires independent
review.
