# Failed / Controlled Negative Tests Archive

This directory stores documented negative test results, failed verification experiments, and boundary breakdown cases:

1. **Zero-Radius Filament Singularity Breakdown (ZT-006.2.19)**:
   - Root cause: Ideal filaments of zero radius yield $B \sim 1/r$, producing diverging local energy density $u \sim B^2 \sim 1/r^2$. The volume integral of $1/r^2$ diverges logarithmically, causing Green integrals to oscillate wildly with grid refinement.
   - Resolution: Physical finite-thickness conductor regularization (ZT-006.2.21).

2. **Uncompensated Celestial Time Jump**:
   - Backward time jump $\Delta t = -1.0\text{ s}$ without frame lock causes traveler to materialize in vacuum $29.8\text{ km}$ away from Earth.

3. **High-Frequency MQS Breakdown**:
   - Drives above $10\text{ MHz}$ ($ka > 0.5$) violate quasistatic limits and trigger `FULL-WAVE REQUIRED`.
