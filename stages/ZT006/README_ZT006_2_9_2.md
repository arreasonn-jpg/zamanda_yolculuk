# ZT-006.2.9.2 — Physics-Corrected Gauss Green Quadrature

ZT-006.2.9.1 exposed the same normalization issue seen previously in the
midpoint branch: the Gauss Green integral returned the raw
∫ Tbar/|x-x'| dV term without the required Einstein prefactor.

This version applies exactly one:

    4G/c^4

factor inside the Gauss Green kernel and adds a regression test.

IMPORTANT:
- ZT-006.2.9.1 h values around 1e-2 were raw integrals, not physical metric
  perturbations, and must not be interpreted as h_mu_nu.
- The corrected output should return to the weak-field scale.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_9_1 --mode quick
python -m zt005.run_zt006_2_9_1 --mode full
```

No Earth source or CTC conclusion is made here.
