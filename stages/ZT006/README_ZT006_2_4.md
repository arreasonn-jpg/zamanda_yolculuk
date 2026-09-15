# ZT-006.2.4 — Finite-Wire Regularization / Source-Singularity Diagnostics

ZT-006.2.3 revealed very different field sensitivity across geometries:
G4 is especially unstable, suggesting that the zero-thickness filament source
model is being sampled too close to the source.

ZT-006.2.4 does not hide this problem. It explicitly measures the fraction of
the quadrature volume within 2, 5, 10, and 20 mm of a source filament and
recomputes h_mu_nu after excluding those points.

Important:
- This is a **numerical regularization diagnostic**, not a claim that the
  physical device has vacuum holes around its conductors.
- A future physical conductor model must replace the zero-thickness filament
  before a final source integral is trusted.
- The purpose here is to determine whether the unstable G2/G3/G4 behavior is
  dominated by near-wire singular behavior.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_4 --mode quick
python -m zt005.run_zt006_2_4 --mode full
```

No Earth/CTC conclusion is made here.
