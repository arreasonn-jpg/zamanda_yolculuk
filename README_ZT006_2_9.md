# ZT-006.2.9 — Frozen T_mu_nu + Deterministic Green-Integral Convergence

ZT-006.2.8 showed that field/T_mu_nu convergence at fixed validation points is
reproducible across random seeds. The remaining uncertainty is therefore the
source integration itself.

This stage computes T_mu_nu on each deterministic cylindrical source grid and
then evaluates the Green integral at fixed observation points. Results are
compared against the finest grid in the same run.

The purpose is to separate:
- EM/T_mu_nu discretization;
- Green-kernel quadrature convergence.

Quick:
(6,12,6) -> (8,16,8) -> (12,24,12)

Full:
(8,16,8) -> (12,24,12) -> (16,32,16) -> (20,40,20)

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_9 --mode quick
python -m zt005.run_zt006_2_9 --mode full
```

No Earth source or CTC conclusion is made.
