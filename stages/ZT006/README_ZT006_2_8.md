# ZT-006.2.8 — Robust Cross-Seed / Distance-Aware Convergence

ZT-006.2.7 showed extremely small field/T_mu_nu errors for all geometries when
evaluated on a fixed point set. This stage checks whether that result is
reproducible across independent point sets.

For each geometry:
- multiple random seeds;
- h = 0.02, 0.01 versus h = 0.005 reference;
- distance cutoffs 0, 10, 20 mm;
- robust mean/median/p95/max relative errors.

This specifically tests whether the earlier G4 sensitivity was caused by an
unrepresentative spatial sample or near-wire outliers.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_8 --mode quick
python -m zt005.run_zt006_2_8 --mode full
```
No GR/CTC conclusion is made here.
