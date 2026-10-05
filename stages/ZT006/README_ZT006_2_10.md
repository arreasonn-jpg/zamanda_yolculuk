# ZT-006.2.10 — Finite-Wire Distance-Cutoff Convergence

ZT-006.2.9.1 showed that changing quadrature families did not restore
monotone convergence. This stage tests a different hypothesis:

The ideal zero-thickness filament model may make the source stress-energy too
singular for direct volume quadrature.

For each geometry and quadrature resolution, points closer than a physical
distance cutoff to any filament sample are excluded:

- 5 mm
- 10 mm
- 20 mm
- 30 mm

This is an explicit numerical/finite-wire regularization diagnostic. It is
not yet a physical conductor model.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_10 --mode quick
python -m zt005.run_zt006_2_10 --mode full
```

No Earth source or CTC conclusion is made here.
