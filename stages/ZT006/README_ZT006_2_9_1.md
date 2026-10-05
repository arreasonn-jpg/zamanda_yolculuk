# ZT-006.2.9.1 — Gauss-Legendre / Periodic-Phi Green Quadrature

ZT-006.2.9 revealed that uniform midpoint quadrature does not converge
monotonically, and it also contained a JSON serialization bug. This stage
fixes the serialization and adds a second deterministic quadrature family.

Source integration uses:
- Gauss-Legendre in radius
- Gauss-Legendre in axial coordinate
- uniform periodic phi

The Jacobian r is included in the weights.

The run also reports a matched midpoint calculation for comparison.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_9_1 --mode quick
python -m zt005.run_zt006_2_9_1 --mode full
```

No Earth source or CTC conclusion is made.
