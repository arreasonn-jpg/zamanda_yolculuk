# Benchmark Spacetimes and Exact Solutions

## 1. Minkowski Spacetime
- Metric: $ds^2 = -dt^2 + dx^2 + dy^2 + dz^2$
- Properties: $\Gamma^\mu_{\alpha\beta} = 0$, $R^\rho{}_{\sigma\mu\nu} = 0$, $G_{\mu\nu} = 0$.
- Geodesics: straight worldlines with constant 4-velocity.

## 2. Schwarzschild Spacetime
- Metric in Schwarzschild coordinates $(t, r, \theta, \phi)$:
  $$ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1}dr^2 + r^2(d\theta^2 + \sin^2\theta d\phi^2)$$
- Properties: Vacuum solution ($R_{\mu\nu} = 0$, $G_{\mu\nu} = 0$), Kretschmann invariant $K = 48 M^2 / r^6$.
- Geodesics: stable circular equatorial orbits exist for $r > 3M$. Validates perihelion precession and RK4 orbital integration.

## 3. Kerr Spacetime
- Metric in Boyer-Lindquist coordinates for spinning black hole ($M, a$).
- Properties: Vacuum solution ($R_{\mu\nu} = 0$, $G_{\mu\nu} = 0$). Frame-dragging effect $g_{t\phi} \ne 0$.
- Geodesics: tests equatorial timelike orbits with non-zero frame-dragging.

## 4. Gödel Universe
- Exact dust solution of Einstein equations with cosmological constant:
  $$ds^2 = \frac{2}{\omega^2} \left[ -dt^2 + dr^2 - (\sinh^4 r - \sinh^2 r) d\phi^2 + dz^2 - 2\sqrt{2}\sinh^2 r dt d\phi \right]$$
- Properties: Homogeneous, rotating dust with energy density $\rho = \omega^2 / (4\pi G)$.
- Closed Timelike Curves (CTCs): $g_{\phi\phi} < 0$ occurs strictly for $r > \text{asinh}(1) \approx 0.8814$.

## 5. Idealized Tipler Cylinder
- Interior metric of an infinitely long, rigidly rotating cylinder of dust:
  $$ds^2 = -dt^2 + dr^2 + r^2(1 - a^2 r^2) d\phi^2 + dz^2 - 2 a r^2 dt d\phi$$
- Properties: Circular azimuthal curves ($t=\text{const}, r=\text{const}, z=\text{const}$) become timelike ($g_{\phi\phi} < 0$) for $a \cdot r > 1$.
