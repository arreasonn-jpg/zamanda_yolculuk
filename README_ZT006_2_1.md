# ZT-006.2.1 — Deterministic Cylindrical Quadrature

The Monte-Carlo source integral in ZT-006.2 did not converge reliably.

This stage replaces it with a deterministic midpoint quadrature in cylindrical
coordinates:

    dV = r dr dphi dz

The Green-function source integral is evaluated as:

    hbar_mu_nu(x) ≈ (4G/c^4) ∑ Tbar_mu_nu(x') w / |x-x'|

using three resolutions:

    (Nr,Nphi,Nz) = (8,16,8)
                         (12,24,12)
                         (16,32,16)

The purpose is to determine whether h00 and h0i stabilize under systematic
resolution refinement.

Run:

```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_1
```

No Earth source, full GR solver, or CTC conclusion is made here.
