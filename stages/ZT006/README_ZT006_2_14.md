# ZT-006.2.14 — Finite-Wire Cross-Section Convergence

Tests convergence of the finite-radius conductor approximation with increasing
cross-section subfilament count.

Quick:
- radii 1, 2.5, 5 mm
- counts 2, 4, 6, 8

Full:
- radii 1, 2.5, 5, 7.5 mm
- counts 2, 4, 6, 8, 12

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_14 --mode quick
python -m zt005.run_zt006_2_14 --mode full
```
No Earth/CTC conclusion is made here.
