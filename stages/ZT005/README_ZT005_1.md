
# ZT-005.1 — Analytic Benchmark & Convergence

ZT-005.1 does not add more physical assumptions. It asks whether the
numerical electromagnetic solver converges toward a known solution.

## Gold-standard benchmark

For a circular filament loop of radius R carrying current I, the magnetic
field on the symmetry axis is

For the code's +t loop orientation, B_y(z) = -mu0 I R^2 / [2 (R^2 + z^2)^(3/2)].

The benchmark checks:

1. Biot-Savart filament discretization convergence.
2. curl(A) spatial-grid convergence against the analytic axis field.

This separates two error sources that were mixed together in ZT-005.

## Acceptance philosophy

No arbitrary "pass" threshold is hard-coded yet. We first inspect the
convergence curve. A reliable solver should show decreasing error as:
- filament segments increase;
- spatial grid resolution increases.

Only after this is demonstrated should G1-G4 Maxwell validation be resumed.

Run:

```powershell
pip install -e .
pytest
python -m zt005.run_zt005_1
```

No CTC or GR conclusion is made here.
