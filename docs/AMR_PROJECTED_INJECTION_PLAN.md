# Plan: constraint projection onto the injected AMR band

Status: FULLY IMPLEMENTED (v13.2.0-v13.4.0). v13.2.0-alpha1:
flat-Laplacian band solve, POSITIVE result (H_f/H0 @ 1000 = 0.9319
vs baseline 1.0756). v13.3.0: exact gamma-tilde covariant operator +
iso-reference metric, NEUTRAL — the gain is operator-robust.
v13.4.0: the momentum Vector-Poisson increment below is also
NEUTRAL on H (band A-divergence cleanable at ratio ~0.43/step but
the Hamiltonian headline is unchanged: 0.9371). Both increments are
measured and closed; see limitations.md v13.2-v13.4 sections.

## Diagnosis chain so far (why projection is the next lever)

- v12.9-v12.15: the fixed-patch AMR deficit is confined to the
  coarse-to-fine boundary injection; neither cubic prolongation
  (v12.11), finer fine resolution (v12.13) nor Gamma-consistent band
  injection (v12.15) moves the Hamiltonian constraint (H_f/H0 at
  step 1000 stays 1.0756).
- v13.0.0: the prolonged fields' second-derivative (Ricci) mismatch
  is the measured carrier; a moving window does not change it.
- Therefore: leave the interpolation machinery untouched and instead
  correct the injected data itself so that it satisfies the fine
  grid's discrete Hamiltonian constraint in the band.

## Mechanism (band-local elliptic conformal solve)

Keep (gtilde, K, Atilde, alpha, beta, Gamma) exactly as prolongated.
Keep phi's boundary values. In the INTERIOR of the injected band,
replace prolongated phi with the solution of the discrete
Hamiltonian equation for the conformal factor:

    Lap(e^{4 phi} - 1) = -(1/8) * <H residual source>

implemented as a bounded number (10-30) of damped Jacobi sweeps on
the band, with Dirichlet boundary values taken from the prolongated
data. The tilde-divergence of Atilde enters the source side of the
equation (the standard conformal-solve route used in puncture and
initial-data codes; here only the band-local version is attempted).

Optional second increment (only if the first shows motion): a local
Vector Poisson correction for Atilde — solve

    Lap(v) = -grad(div Atilde)

on the band and subtract the resulting gradient-sourced v field from
Atilde, again by Jacobi sweeps.

## Experiment protocol (identical to the v12.x harness)

- The same binary, grids, dt, injection layers (N_c = 17, N_f = 17,
  factor = 2, dt = 2e-3, 4 layers, 1000 steps).
- Baseline: the current injection (H_f/H0 at step 1000 = 1.0756).
- Success criteria (any of): material reduction of H_f/H0 at step
  1000, materially slower growth between steps 200-1000, or a
  smaller H_f growth rate with the band-solve enabled.
- Failure is also a result: recorded in limitations.md like its
  predecessors.

## Cost estimate

Band-local only: the fine window is 17^3 and the Jacobi sweeps run
on the band subset; this is cheap numpy slicing work and does not
touch the RK4/RHS numerics. No baseline drift is expected, since
only injected band values are modified.

## Related performance note (separate; do not mix in)

Profiling the 4x AMR step (2026-10-20) showed that `numpy.einsum`
calls inside `bssn_rhs_vec.py` dominate the fine-grid step
(10.8 s of a 20.9 s profile at N=33). Adding `optimize=True` would
likely speed up fine-grid evolution 3-5x, but shifts the summation
order (round-off-level numerical changes). It must NOT be switched
inside a result run; schedule as its own version with a full-suite
re-validation and the baseline reproduction check (the 1.0756
reference must still reproduce within tolerance, or the change is
documented as a version-tagged numerical change).
