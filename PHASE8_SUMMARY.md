# Phase 8: Energy Budget Skeleton - COMPLETED

## Summary
Successfully implemented the energy budget skeleton (Phase 8 priority from roadmap) that translates the h ~ 10⁻⁴⁶ metric perturbation scale into explicit energy requirements.

## Files Created
- `src/zt005/energy_budget.py` - Energy budget computation module
- `results/validated/energy_budget.json` - Computed energy budgets for all configs

## Results

| Config | P_in (W) | P_loss (W) | P_stored (J) | η | h_current | Gap (orders) |
|--------|----------|------------|--------------|-----|-----------|--------------|
| G1 | 84.0 | 84.0 | 7.37e-06 | 8.78e-08 | 6.09e-50 | 46.2 |
| G2 | 84.0 | 84.0 | 3.66e-06 | 4.35e-08 | 3.02e-50 | 46.5 |
| G3 | 84.0 | 84.0 | 1.85e-08 | 2.20e-10 | 1.53e-52 | 48.8 |
| G4 | 84.0 | 84.0 | 1.04e-06 | 1.24e-08 | 8.61e-51 | 47.1 |

## Key Findings
1. **Input power**: ~84 W for all configs (100 A peak, 1 MHz, 1 mm² copper wire)
2. **Efficiency**: Extremely low (10⁻⁸ to 10⁻¹⁰) - almost all power goes to ohmic losses
3. **Stored EM energy**: Microjoules (10⁻⁶ to 10⁻⁸ J)
4. **Metric perturbation**: h ~ 10⁻⁵⁰ to 10⁻⁵² (far below CTC-relevant h ~ 10⁻³)
5. **Gap**: 46-49 orders of magnitude to reach h = 10⁻³

## MQS Limitation Documented
Created `docs/mqs_limitation.md` documenting the fundamental limitation:
- MQS approximation neglects displacement current
- Stress-energy tensor is NOT conserved (∂_μ T^μν ≠ 0)
- Lorenz gauge residual fails for device source (scaled_residual ~ 1.0)
- Manufactured conserved source passes (scaled_residual ~ 0.03)
- This is a known limitation, not a bug

## Test Status
✅ All 779 tests pass (12 min runtime)

## Next Priorities (from roadmap)
1. **BSSN AMR projection long runs (v14.1)** - Critical path for inspiral marathon
2. **BSSN ↔ retarded_hbar cross-check** - Numeric cross-validation
3. **Full-wave Maxwell (ZT-007)** - Deferred until frequency scans demand it