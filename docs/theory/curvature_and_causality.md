# Spacetime Curvature, Geodesics, and Causal Structure

## 1. Riemann Curvature Chain
From any spacetime metric $g_{\mu\nu}$, the curvature is derived through the standard differential geometric sequence:

1. Metric inverse: $g^{\mu\nu} = (g^{-1})_{\mu\nu}$
2. Christoffel symbols (Levi-Civita connection):
   $$\Gamma^\mu_{\alpha\beta} = \frac{1}{2} g^{\mu\lambda} \left( \partial_\alpha g_{\beta\lambda} + \partial_\beta g_{\alpha\lambda} - \partial_\lambda g_{\alpha\beta} \right)$$
3. Riemann curvature tensor:
   $$R^\rho{}_{\sigma\mu\nu} = \partial_\mu \Gamma^\rho_{\nu\sigma} - \partial_\nu \Gamma^\rho_{\mu\sigma} + \Gamma^\rho_{\mu\lambda} \Gamma^\lambda_{\nu\sigma} - \Gamma^\rho_{\nu\lambda} \Gamma^\lambda_{\mu\sigma}$$
4. Ricci tensor:
   $$R_{\mu\nu} = R^\lambda{}_{\mu\lambda\nu}$$
5. Ricci scalar:
   $$R = g^{\mu\nu} R_{\mu\nu}$$
6. Einstein tensor:
   $$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} g_{\mu\nu} R$$

## 2. Einstein Field Equations and Residual
$$G_{\mu\nu} = \frac{8\pi G}{c^4} T_{\mu\nu}$$
In the linearized perturbation regime $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$, the solver is strictly classified as a **"linearized metric perturbation solver"** rather than a full nonlinear GR solver until self-gravitating backreaction is solved iteratively.

## 3. Geodesic Motion and Proper Time
The worldline of a test body is determined by the geodesic equation:
$$\frac{d^2 x^\mu}{d\tau^2} + \Gamma^\mu_{\alpha\beta} \frac{dx^\alpha}{d\tau} \frac{dx^\beta}{d\tau} = 0$$
The metric norm $g_{\mu\nu} u^\mu u^\nu = -c^2$ is an exact constant of motion along timelike geodesics.

## 4. Proper Time vs. Coordinate Time
The traveler's internal clock measures proper time $\tau$:
$$c^2 d\tau^2 = - g_{\mu\nu} dx^\mu dx^\nu$$
while external observers measure coordinate time $t = x^0 / c$.
Backward time travel requires:
$$\Delta t = t_{\rm final} - t_{\rm initial} < 0 \quad \text{while} \quad \Delta \tau > 0$$
Along closed timelike curves (CTCs), a timelike curve returns to its starting spacetime coordinate:
$$x^\mu(\lambda_1) = x^\mu(\lambda_0) \quad \text{with} \quad ds^2 < 0 \text{ everywhere.}$$
