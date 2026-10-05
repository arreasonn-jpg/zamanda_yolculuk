"""Manufactured stress-energy sources for conservation/gauge validation.

The key validation problem is: we need a stress-energy distribution that is
EXACTLY conserved (partial_mu T^{mu nu} = 0 pointwise), smooth, and
compactly supported, so that numerical residuals measure the solver and not
the source.

Construction used here
----------------------
Take any smooth bump phi(x) and define

    T^{mu nu} := G^{(1) mu nu}[h],      h_munu = phi * eta_munu,

the linearized Einstein tensor of a pure-trace perturbation.  By the
linearized Bianchi identity, partial_mu G^{(1) mu nu} = 0 holds IDENTICALLY
for any phi, so this source is conserved to machine precision before any
discretization.  For h_munu = phi eta_munu the tensor simplifies to

    T^{mu nu} = - partial^mu partial^nu phi + eta^{mu nu} Box phi

(all derivatives analytically known for a Gaussian phi).

A deliberately BROKEN source (T^{00} = phi, all other components zero) is
also provided; any honest conservation diagnostic must flag it.

Coordinates x^mu = (c t, x, y, z), signature (-, +, +, +), SI units.
"""

from __future__ import annotations

import numpy as np

C = 299_792_458.0


class GaussianBump:
    """phi(t, x) = A exp(-Q), Q = (t-t0)^2/(2 st^2) + |x-x0|^2/(2 sx^2)."""

    def __init__(self, A=1.0, t0=0.0, x0=(0.0, 0.0, 0.0),
                 sigma_t=1e-7, sigma_x=1.0):
        self.A = float(A)
        self.t0 = float(t0)
        self.x0 = np.asarray(x0, float)
        self.st = float(sigma_t)
        self.sx = float(sigma_x)

    def value(self, t, x):
        t = np.asarray(t, float)
        x = np.asarray(x, float)
        dt = t - self.t0
        dx = x - self.x0
        Q = dt * dt / (2 * self.st ** 2) + np.sum(dx * dx, axis=-1) \
            / (2 * self.sx ** 2)
        return self.A * np.exp(-Q)

    # --- analytic partial derivatives (coordinates (ct, x, y, z)) -----
    def d0(self, t, x):
        """partial_0 phi = (1/c) d phi / dt"""
        dt = np.asarray(t, float) - self.t0
        return -self.value(t, x) * dt / self.st ** 2 / C

    def d00(self, t, x):
        phi = self.value(t, x)
        dt = np.asarray(t, float) - self.t0
        return phi * ((dt / self.st ** 2) ** 2 - 1.0 / self.st ** 2) / C ** 2

    def di(self, t, x):
        """spatial gradient, shape (..., 3)"""
        phi = self.value(t, x)
        dx = np.asarray(x, float) - self.x0
        return -phi[..., None] * dx / self.sx ** 2

    def d0i(self, t, x):
        """partial_0 partial_i phi = (1/c) d_t d_i phi, shape (..., 3)"""
        phi = self.value(t, x)
        dt = np.asarray(t, float) - self.t0
        dx = np.asarray(x, float) - self.x0
        return (phi[..., None] * dt[..., None] * dx
                / (self.st ** 2 * self.sx ** 2) / C)

    def dij(self, t, x):
        """spatial Hessian, shape (..., 3, 3)"""
        phi = self.value(t, x)
        dx = np.asarray(x, float) - self.x0
        eye = np.eye(3)
        return (phi[..., None, None]
                * (dx[..., :, None] * dx[..., None, :] / self.sx ** 4
                   - eye / self.sx ** 2))

    def laplacian(self, t, x):
        phi = self.value(t, x)
        dx = np.asarray(x, float) - self.x0
        r2 = np.sum(dx * dx, axis=-1)
        return phi * (r2 / self.sx ** 4 - 3.0 / self.sx ** 2)

    def dalembertian(self, t, x):
        """Box phi = -partial_0^2 phi + laplacian phi   [(-+++)]"""
        return -self.d00(t, x) + self.laplacian(t, x)


def conserved_T(bump, t, x):
    """T^{mu nu} = -partial^mu partial^nu phi + eta^{mu nu} Box phi.

    Vectorized over source points; ``t`` may be per-point (shape (Ns,)).
    Returns (Ns, 4, 4).  All derivatives are computed from a single
    evaluation of phi for speed.
    """
    t = np.atleast_1d(np.asarray(t, float))
    x = np.asarray(x, float)
    if x.ndim == 1:
        x = np.broadcast_to(x, t.shape + (3,)).copy()
    n = len(t)

    dt_ = (t - bump.t0)[..., None]
    dx = x - bump.x0
    r2 = np.sum(dx * dx, axis=-1)
    Q = dt_[..., 0] ** 2 / (2 * bump.st ** 2) + r2 / (2 * bump.sx ** 2)
    phi = bump.A * np.exp(-Q)

    st2, sx2 = bump.st ** 2, bump.sx ** 2
    d00 = phi * ((dt_[..., 0] / st2) ** 2 - 1.0 / st2) / C ** 2
    d0i = phi[..., None] * dt_ * dx / (st2 * sx2) / C
    dii = (phi[..., None, None]
           * (dx[..., :, None] * dx[..., None, :] / sx2 ** 2
              - np.eye(3) / sx2))
    lap = phi * (r2 / sx2 ** 2 - 3.0 / sx2)
    box = -d00 + lap

    T = np.zeros((n, 4, 4))
    # raising/lowering with eta: partial^0 = -partial_0, partial^i = partial_i
    T[:, 0, 0] = -d00 - box          # -(+d00) + eta^00 box, eta^00=-1
    T[:, 0, 1:] = d0i                # -partial^0 partial^i = +partial_0 partial_i
    T[:, 1:, 0] = d0i
    T[:, 1:, 1:] = -dii + box[:, None, None] * np.eye(3)
    return T


def conserved_source_fn(bump):
    """source_fn(t_array, x_points) adapter for retarded_gr.retarded_hbar."""
    def fn(t, x):
        return conserved_T(bump, t, x)
    return fn


def broken_T(bump, t, x):
    """NOT conserved control source: T^{00} = phi and isotropic spatial
    stress T^{ij} = delta^{ij} phi, but no momentum density T^{0i}.

    Then partial_i T^{ij} = partial_j phi =/= 0, a violation at the
    spatial-gradient scale of the bump (a strong, unambiguous signal).
    """
    t = np.atleast_1d(np.asarray(t, float))
    x = np.asarray(x, float)
    if x.ndim == 1:
        x = np.broadcast_to(x, t.shape + (3,)).copy()
    phi = bump.value(t, x)
    T = np.zeros((len(t), 4, 4))
    T[:, 0, 0] = phi
    T[:, 1, 1] = phi
    T[:, 2, 2] = phi
    T[:, 3, 3] = phi
    return T


def broken_source_fn(bump):
    def fn(t, x):
        return broken_T(bump, t, x)
    return fn
