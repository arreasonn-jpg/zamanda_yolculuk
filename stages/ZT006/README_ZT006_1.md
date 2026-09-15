# ZT-006.1 — Spatial Linearized-GR Green-Function Diagnostic

ZT-006.1 replaces the single-distance characteristic estimate with a spatial
volume integral over the active cylinder:

    hbar_mu_nu(x) ≈ (4G/c^4) ∫ [Tbar_mu_nu(x') / |x-x'|] d^3x'

The source is sampled uniformly in the active cylinder using Monte Carlo
points. Eight observation points are evaluated inside the active region.

A small kernel softening (0.5h) is used only to regularize the discrete
sampling of the 1/r Green function. This is a numerical regularization, not
a physical new interaction.

After obtaining hbar, the physical perturbation is reconstructed via

    h_mu_nu = hbar_mu_nu - 1/2 eta_mu_nu hbar

The calculation is quasi-static and linearized. It is NOT a full numerical
solution of the Einstein equations, and it does not test CTC existence.

Run:

```powershell
pip install -e .
pytest
python -m zt005.run_zt006_1
```
