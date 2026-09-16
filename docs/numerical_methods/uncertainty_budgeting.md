# Uncertainty Budgeting Formulation

Every computed metric perturbation component must be reported with explicit uncertainty bounds:
$$h = h_{\rm estimate} \pm \delta h_{\rm total}$$

## Six Independent Error Sources
1. **$\delta_{\rm grid}$ (Grid Discretization Error)**:
   Estimated from step difference between finest grid levels $N=36$ and $N=28$:
   $$\delta_{\rm grid} = |h(36) - h(28)|$$
2. **$\delta_{\rm quad}$ (Quadrature Scheme Error)**:
   Discrepancy between Solver A (Gauss-cylindrical) and Solver B (adaptive Cartesian):
   $$\delta_{\rm quad} = |h_A - h_B|$$
3. **$\delta_{\rm src}$ (Source Segment Discretization Error)**:
   Discrepancy between filament circle/helix segment resolutions ($n=128$ vs $n=256$).
4. **$\delta_{\rm wire}$ (Finite Conductor Discretization Error)**:
   Variation from cross-sectional wire subdivision ($(n_r, n_\phi) = (2, 8)$ vs $(4, 12)$).
5. **$\delta_{\rm fp}$ (Floating Point Roundoff Error)**:
   IEEE-754 double precision accumulation:
   $$\delta_{\rm fp} \approx \epsilon_{\rm mach} \cdot \kappa \cdot |h| \sim 10^{-14} |h|$$
6. **$\delta_{\rm model}$ (Physical Model Approximation Error)**:
   Trunctation of MQS $\mathcal{O}((ka)^2)$ and linearized GR backreaction $\mathcal{O}(h^2)$:
   $$\delta_{\rm model} \approx (ka)^2 |h| + |h|^2$$

## Combined Uncertainty
The total uncertainty is calculated via root-sum-square (RSS) quadrature:
$$\delta h_{\rm total} = \sqrt{\delta_{\rm grid}^2 + \delta_{\rm quad}^2 + \delta_{\rm src}^2 + \delta_{\rm wire}^2 + \delta_{\rm fp}^2 + \delta_{\rm model}^2}$$
Typical baseline result:
$$h_{00} = (1.980000 \pm 0.008388) \times 10^{-46} \quad [\pm 0.42\%]$$
