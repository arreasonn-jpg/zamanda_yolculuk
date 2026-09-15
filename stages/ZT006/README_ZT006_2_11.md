# ZT-006.2.11 — Finite-Cross-Section Conductor Model

ZT-006.2.10 showed that the minimum distance from quadrature nodes to the
ideal filaments changes strongly with grid refinement. This stage replaces
the zero-thickness filament with a finite-radius approximation: the same total
current is distributed over multiple parallel offset filament paths across a
circular conductor cross-section.

Quick radii: 0, 1, 2.5, 5 mm.
Full adds: 7.5 mm.

The goal is to determine whether field/stress-energy behavior stabilizes when
the source singularity is physically regularized.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_11 --mode quick
python -m zt005.run_zt006_2_11 --mode full
```

This is a numerical finite-wire approximation, not yet a material/thermal/
eddy-current model. No Earth/CTC conclusion is made here.
