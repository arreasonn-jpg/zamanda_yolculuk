"""
ZT-006.2.21 — Axisymmetric Poisson Solver for Linearized Einstein Eq.

Solves:
    ∇²h = -(16πG/c⁴) T00

in cylindrical coordinates (r, z), axisymmetric.

This replaces the Green-function integral approach used in ZT-006.2.x,
which was numerically unstable for G2/G3/G4 geometries.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

# Physical constants (SI)
G = 6.67430e-11
C = 2.99792458e8
EINSTEIN_FACTOR = 16.0 * np.pi * G / C**4   # ≈ 3.305e-43


def solve_poisson_axisym(T00, r_vals, z_vals):
    """
    Solve axisymmetric Poisson equation in (r, z):
        ∂²h/∂r² + (1/r) ∂h/∂r + ∂²h/∂z² = -S(r, z)
    where S = EINSTEIN_FACTOR * T00.

    Boundary conditions:
        r = 0   : Neumann (even symmetry, ∂h/∂r = 0)
        r = R   : Dirichlet h = 0
        z = ±Z  : Dirichlet h = 0

    Parameters
    ----------
    T00 : (Nr, Nz) ndarray
    r_vals : (Nr,) ndarray
    z_vals : (Nz,) ndarray

    Returns
    -------
    h : (Nr, Nz) ndarray
    """
    Nr, Nz = T00.shape
    dr = r_vals[1] - r_vals[0]
    dz = z_vals[1] - z_vals[0]
    inv_dr2 = 1.0 / dr**2
    inv_dz2 = 1.0 / dz**2

    N = Nr * Nz
    rows, cols, vals, b = [], [], [], np.zeros(N)

    for i in range(Nr):
        for j in range(Nz):
            k = i * Nz + j

            # --- Dirichlet at outer boundaries ---
            if i == Nr - 1 or j == 0 or j == Nz - 1:
                rows.append(k); cols.append(k); vals.append(1.0)
                continue

            if i == 0:
                # Axis (r=0): use L'Hôpital limit
                # (1/r) ∂h/∂r  →  ∂²h/∂r²
                # So equation becomes: 2∂²h/∂r² + ∂²h/∂z² = -S
                # FD with even symmetry h[-1,j] = h[1,j]:
                rows.append(k); cols.append(k);                vals.append(-4.0*inv_dr2 - 2.0*inv_dz2)
                rows.append(k); cols.append((i+1)*Nz + j);      vals.append( 4.0*inv_dr2)
                rows.append(k); cols.append(i*Nz + (j+1));      vals.append( inv_dz2)
                rows.append(k); cols.append(i*Nz + (j-1));      vals.append( inv_dz2)
            else:
                r = r_vals[i]
                coeff_rp = inv_dr2 + 1.0 / (2.0 * r * dr)
                coeff_rm = inv_dr2 - 1.0 / (2.0 * r * dr)

                rows.append(k); cols.append(k);                vals.append(-2.0*inv_dr2 - 2.0*inv_dz2)
                rows.append(k); cols.append((i+1)*Nz + j);      vals.append(coeff_rp)
                rows.append(k); cols.append((i-1)*Nz + j);      vals.append(coeff_rm)
                rows.append(k); cols.append(i*Nz + (j+1));      vals.append(inv_dz2)
                rows.append(k); cols.append(i*Nz + (j-1));      vals.append(inv_dz2)

            b[k] = -EINSTEIN_FACTOR * T00[i, j]

    A = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))
    h_flat = spsolve(A, b)
    return h_flat.reshape(Nr, Nz)