# ZT-006.2.3 — T_mu_nu / Green-Kernel / Softening Error Budget

This stage separates three numerical sensitivities:

1. **Field reconstruction sensitivity**: T_mu_nu at fixed spatial points as
   the 4th-order curl step h changes.
2. **Green-kernel resolution sensitivity**: the spatial source integration at
   increasing cylindrical quadrature resolutions.
3. **Kernel softening sensitivity**: h_mu_nu response to 0, 1, 2.5, 5, 10 mm
   regularization of the 1/r kernel.

This is an error-budget study. It does not treat any one component as fully
validated merely because another component converges.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_3 --mode quick
python -m zt005.run_zt006_2_3 --mode full
```

No Earth source, exact curved-spacetime solver, or CTC conclusion is made here.
