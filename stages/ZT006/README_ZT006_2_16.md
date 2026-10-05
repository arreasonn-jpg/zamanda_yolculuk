# ZT-006.2.16 — Source-Model Gate + Linearized Metric Validation

ZT-006.2.15 establishes stable area-weighted finite-wire cross-section convergence.
This stage begins the transition from source validation to metric validation.

It evaluates h_mu_nu at several observation points, constructs g=eta+h, checks
metric determinant/eigenvalue sanity, and checks a rest-frame timelike vector.
This is not a CTC detector and makes no time-travel claim.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_16 --mode quick
python -m zt005.run_zt006_2_16 --mode full
```
