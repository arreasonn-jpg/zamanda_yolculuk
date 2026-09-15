# Numerical methods

## Quadrature of 1/R kernels

Midpoint (cell-collocation) quadrature over structured grids:

    ∫ f(x′)/|x−x′| d³x′ ≈ Σ_j w_j f(x_j)/|x−x_j|

Known defect and its fix (ZT-006.3): when an observer sits on a source
point, w/R diverges while the true cell integral is finite. The solver
replaces coincident cells by the spherical-average kernel
⟨1/R⟩ = 3/(2 r_cell), r_cell = (3w/4π)^{1/3}
(`retarded_gr.singular_cell_kernel`). Without this correction a single
coincident cell produces a spurious spike ~10¹²× the physical value and
destroys every downstream derivative.

Convergence status of the *frozen* kernel (ZT-006.2.19/.20): **NOT
converged**, see results/archived and stages/ZT006/README.md. The
retarded solver inherits the same quadrature family; its convergence
discipline must be re-established per configuration before any physics
claim.

## Finite differences

Second-order central differences on regular spacetime grids for

* ∂_μ T^{μν} (stress_energy_conservation.divergence_residual)
* ∂_μ h̄^{μν} (gauge_checks.lorenz_residual)
* Christoffel symbols from metric callables
  (benchmark_metrics.Metric.christoffel_numeric)

Critical implementation detail: the divergence along coordinate μ acts on
the **matching tensor row** (∂_x differentiates T^{1ν}, not T^{μν}
componentwise). Getting the row index wrong produces large structured
residuals even for exactly conserved sources — this bug was found and
fixed in ZT-006.3.

Error hierarchy that governs residual tolerances:

    discretization error of residual
        ~ O(h_obs²) · third derivatives of the field
        + O(h_src²) · quadrature noise amplified by 1/h_obs

Hence residual tests must state grid spacings together with tolerances.

## Time integration

Fixed-step RK4 for geodesics:

    d²x^μ/dτ² = −Γ^μ_{αβ} ẋ^α ẋ^β

Quality diagnostic: drift of the exact invariant g_μν u^μ u^ν
(< 1e−9 with analytic connections, < 1e−5 with numeric connections in
the current test suite).

## Manufactured-source method

Validation sources are built as

    T^{μν} := G^{(1)μν}[φ η]  =  −∂^μ∂^ν φ + η^{μν} □ φ

with φ a Gaussian bump. By the linearized Bianchi identity this source is
divergence-free *identically*, smooth, and effectively compactly
supported — so measured residuals belong to the solver, not the source.
A deliberately broken control source (isotropic stress without momentum
density) checks that the diagnostic actually fires.

## Known-answer tests

CTC search is validated against exact thresholds (Gödel r > asinh 1,
Tipler ar > 1) and must return *no* CTC for Minkowski / Schwarzschild
exterior / Morris–Thorne. Schwarzschild perihelion advance is compared
against the 1PN prediction 6πM/(a(1−e²)) with measured (a, e).
