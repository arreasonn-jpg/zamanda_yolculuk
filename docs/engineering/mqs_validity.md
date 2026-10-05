# Magnetoquasistatic (MQS) Validity Criteria

The dimensionless parameter governing the validity of the magnetoquasistatic approximation is:
$$ka = \frac{\omega a}{c} = \frac{2\pi f a}{c} = 2\pi \frac{a}{\lambda}$$
where $f$ is the drive frequency, $a$ is the characteristic dimension of the apparatus (e.g. radius $1.0\text{ m}$), and $\lambda = c/f$ is the free-space electromagnetic wavelength.

## Regime Boundaries
1. **$ka < 0.1$ $\implies$ MQS VALID**:
   Wave propagation delay across the apparatus is negligible ($t_{\rm transit} / T_{\rm period} < 1.6\%$). Retardation and dipole radiation losses are negligible.
   * At $f = 1.0\text{ MHz}$, $a = 1.0\text{ m}$: $ka \approx 0.0210 < 0.1$ (**MQS VALID**).
2. **$0.1 \le ka \le 0.5$ $\implies$ MQS MARGINAL**:
   Phase shift across the coil assembly is noticeable. Full-wave correction terms should be tracked.
3. **$ka > 0.5$ $\implies$ FULL-WAVE REQUIRED**:
   The apparatus dimension is comparable to or exceeds the wavelength. Radiation damping and retardation cannot be neglected. Full-wave Maxwell solver required (ZT-007).
