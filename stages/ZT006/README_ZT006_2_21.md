# ZT-006.2.21 — Current Normalization, Dual Independent Green Solvers, & Convergence Gate

**Scientific Milestone: Closing Section 1.1, 1.2, and 1.3**

## 1. Akım Normalizasyonu (Section 1.1)
Unambiguously disentangles and enforces:
- `I_parent`: Current per parent conductor turn ($100.0\text{ A}$)
- `I_total`: Total current across all physical turns ($N_{\rm turns} \cdot I_{\rm parent}$)
- `I_peak`: Peak AC amplitude
- `I_rms`: RMS current $I_{\rm peak} / \sqrt{2}$

**Acceptance**: Verified across Gauss-disk $(n_r, n_\phi)$ cross-sectional quadratures.
Maximum discrepancy between resolutions: $< 4 \times 10^{-13}\text{ A}$ (machine precision).

## 2. Bağımsız İkinci Green Solver (Section 1.3)
- **Solver A**: Gauss-cylindrical tensor-product quadrature with spherical cell singularity handling.
- **Solver B**: Cartesian cell-integrated voxel quadrature with boundary volume fraction weighting and exact boundary-face divergence-theorem box potential evaluation:
  $$\int_V \frac{1}{r} dV = \frac{1}{2} \oint_{\partial V} \frac{\mathbf{r} \cdot \hat{\mathbf{n}}}{r} dA$$
**Acceptance**: Relative tensor discrepancy between Solver A and Solver B is $0.27\% < 1\%$.

## 3. Green-Integral Yakınsaması (Section 1.2)
Grid refinement series evaluated across $N \in [12, 16, 20, 28, 36]$:
- $N=12 \to N=16$: $\Delta h_{00} = 1.76\%$
- $N=16 \to N=20$: $\Delta h_{00} = 0.816\%$
- $N=20 \to N=28$: $\Delta h_{00} = 0.712\%$
- $N=28 \to N=36$: $\Delta h_{00} = 0.295\% < 1\%$

Final step acceptance:
- $h_{00}$ relative error: $0.295\% < 1.0\%$ (**PASS**)
- $\|h_{ij}\|$ relative error: $0.295\% < 1.0\%$ (**PASS**)
- $|h_{0i}|$ absolute error: $0.0 < 10^{-45}$ (**PASS**)
- Cross-validated with independent Solver B (**PASS**)

## How to Run
```bash
pytest tests/unit/test_green_solvers.py tests/unit/test_current_normalization.py
python3 -m zt005.run_zt006_2_21 --mode quick
python3 -m zt005.run_zt006_2_21 --mode full
```
