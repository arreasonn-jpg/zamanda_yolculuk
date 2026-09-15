# ZT-005.4 — Validated Pointwise EM → Stress-Energy

ZT-005.3.2 showed that the pointwise 4th-order curl(A) reconstruction has
much smaller error than the 2nd-order reconstruction at h=0.005 m.

ZT-005.4 uses that validated pointwise field representation to calculate:

- electromagnetic energy density
- Poynting vector
- Maxwell stress tensor
- electromagnetic stress-energy tensor T_mu_nu
- Einstein-source scale (8*pi*G/c^4) T_mu_nu
- angular-momentum-density proxy r x (S/c^2)

This is still a source-side calculation. It does not solve the Einstein field
equations and does not claim metric perturbation or CTC formation.

Run:

```powershell
pip install -e .
pytest
python -m zt005.run_zt005_4
```
