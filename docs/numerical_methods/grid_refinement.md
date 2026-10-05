# Grid Refinement Series and Convergence Methodology

## Grid Refinement Protocol
To eliminate false convergence claims, a systematic grid refinement series is evaluated:
$$N \in [12, 16, 20, 28, 36]$$
where the number of cylindrical nodes is $(n_r, n_\phi, n_z) = (N, 2N, N)$, scaling total nodes from $3,456$ up to $93,312$.

## Sensitivity Dimensions
Convergence is analyzed across three independent sensitivity dimensions:
1. **Observer position**: Center $(0, 0, 0)$, interior radial $(0.3, 0, 0)$, interior axial $(0, 0, 0.5)$, edge $(0.8, 0, 0)$, and exterior $(0, 0, 2.0)$.
2. **Quadrature order**: Comparing varied Legendre polynomial degrees and Fourier segment counts.
3. **Regularization scale**: Regularization radius $\rho$ matched adaptively to local cell volume.

## Acceptance Criteria
Between the final two resolutions ($N=28$ vs $N=36$):
- $h_{00}$ relative difference:
  $$\frac{|h_{00}(36) - h_{00}(28)|}{|h_{00}(36)|} < 1\%$$
- $\|h_{ij}\|$ spatial relative Frobenius difference:
  $$\frac{\|h_{ij}(36) - h_{ij}(28)\|}{\|h_{ij}(36)\|} < 1\%$$
- $|h_{0i}|$ cross-coupling error:
  Absolute difference $< 10^{-45}$ or relative difference $< 1\%$.

Observed convergence on benchmark conserved source:
- $N=12 \to 16$: $1.76\%$
- $N=16 \to 20$: $0.816\%$
- $N=20 \to 28$: $0.712\%$
- $N=28 \to 36$: $0.295\% < 1.0\%$ (**PASS**).
