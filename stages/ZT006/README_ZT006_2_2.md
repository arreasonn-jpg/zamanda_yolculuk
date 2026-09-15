
# ZT-006.2.2 — Deterministic Quadrature Physics-Corrected

ZT-006.2.1 exposed a critical bookkeeping error: the deterministic quadrature
returned the raw Green-function integral

    ∫ Tbar_mu_nu / |x-x'| dV'

but the Einstein-field weak-source equation also requires the prefactor

    4G/c^4.

The previous h00 values (~1e-2) were therefore **not metric perturbations**.
They were unscaled source integrals.

ZT-006.2.2 applies the missing factor 4G/c^4 inside the quadrature routine and
adds a regression test for the prefactor.

This is a mandatory correction before any physical interpretation of h_mu_nu.

Run:

```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_1
```

Expected physical scale should return to the extremely small weak-field regime
seen in ZT-006.1, subject to resolution differences.

No Earth model or CTC conclusion is made here.
