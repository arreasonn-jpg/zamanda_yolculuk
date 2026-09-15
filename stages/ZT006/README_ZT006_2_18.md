ZT-006.2.18 — singular-cell regularized Green-kernel validation.
Because the observer lies inside the source volume, 1/r is integrable but poorly resolved by naive sampling. A near spherical cell uses the exact identity integral(dV/r)=2*pi*rho^2 with local weighted T, while the far field is quadrature-integrated.
Run: python -m pytest; python -m zt005.run_zt006_2_18 --mode quick; python -m zt005.run_zt006_2_18 --mode full
