# MQS Approximation Limitation: Non-Conserved Stress-Energy Tensor

## Summary
The magneto-quasi-static (MQS) approximation used in the device EM stress-energy tensor computation **does not satisfy the full Maxwell equations**, leading to a non-conserved stress-energy tensor (∂_μ T^μν ≠ 0). This causes the Lorenz gauge condition (∂_μ h̄^μν = 0) to fail in the retarded linearized Einstein solver.

## Root Cause
The MQS approximation:
- **Includes**: B = ∇ × A, E = -∂A/∂t (Faraday's law)
- **Neglects**: Displacement current μ₀ε₀∂E/∂t in Ampere's law

This means the EM fields don't satisfy the full Maxwell equations, so the stress-energy tensor T^μν computed from these fields is not conserved.

## Evidence
| Source Type | Scaled Residual | Status |
|-------------|----------------|--------|
| Device (MQS) | 0.82 - 1.13 | **FAIL** |
| Manufactured (Conserved) | 0.03 | **PASS** |

The manufactured source (GaussianBump with T^μν = -∂^μ∂^νφ + η^μν□φ) is analytically conserved and passes the gauge check with scaled_residual ≈ 0.03.

## Impact
- The retarded linearized Einstein solver assumes a conserved source for Lorenz gauge consistency
- Device pipeline produces h̄^μν that doesn't satisfy ∂_μ h̄^μν = 0
- This is a **fundamental limitation of the MQS approximation**, not a bug in the solver

## Options Forward
1. **Accept limitation** (current): Document as known issue, proceed with other priorities
2. **Full-wave Maxwell (ZT-007)**: Implement full Maxwell solver with displacement current
3. **Correction terms**: Add post-hoc correction to enforce conservation

## Roadmap Reference
- ZT-007 (Full-wave Maxwell) is deferred until frequency scans demand it
- Current Phase 8 priority: Energy bookkeeping skeleton
- This limitation should be documented in validation reports

## Files
- `scripts/gauge_audit_device.py` - Audit script
- `src/zt005/batch_em_source.py` - Device EM source (MQS)
- `src/zt005/validated_stress_energy.py` - Stress-energy computation
- `src/zt005/manufactured_sources.py` - Conserved source for validation