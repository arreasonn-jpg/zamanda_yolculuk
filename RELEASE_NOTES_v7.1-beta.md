# v7.1-beta Release Notes

**Date**: 2026-10-31  
**Tag**: `v7.1-beta`  
**Commit**: `v7.1-beta consolidation`

---

## Summary

v7.1-beta completes the refined Kerr + complex scalar + PTUM validation with NEC violation stable on 2:1 AMR grid, and documents the Ford-Roman quantum inequality violation at strong NEC. The two-way coupling with peak-preserving restrict was tested and found to lose NEC violation (Phi decays to zero), confirming one-way coupling as the canonical approach.

---

## Validation Results

### 1. Refined Kerr + Complex + PTUM (100 steps, A_fuel=3.0, a=0.9)

| Metric | Value | Status |
|--------|-------|--------|
| **NEC violation** | ρ_min = -1.5658e+01 (persistent) | ✅ |
| **Negative cells** | 27 (constant) | ✅ |
| **K_max (coarse)** | 0.2485 → 0.2419 (decreasing) | ✅ |
| **K_max (fine)** | 0.3646 → 0.3603 (decreasing) | ✅ |
| **Phi_max** | 3.0002 → 3.0008 (stable) | ✅ |
| **psi_max** | 2.8284 (constant) | ✅ |
| **alpha_min** | 0.8432 (constant, > 0.8) | ✅ |
| **All fields finite** | Yes | ✅ |

**Grid**: Coarse N=9, L=20; Fine N=9, L=10; offset=(2,2,2)  
**Time step**: dt=2e-4, 100 steps  
**Runtime**: ~222s

### 2. Ford-Roman Quantum Inequality Check

**Configuration**: Same as above, N_field_components=2 (PTUM has two real scalars: Phi, psi)

| τ | Massless Ratio | Massless Violates | Massive (μ=1) Ratio | Massive Violates |
|---|----------------|-------------------|---------------------|------------------|
| 0.01 | 0.0000 | False | 0.0000 | False |
| 0.05 | 0.0005 | False | 0.0005 | False |
| 0.10 | 0.0039 | False | 0.0039 | False |
| 0.50 | 0.4913 | False | 0.7677 | False |
| 1.00 | **3.9325** | **True** | **15.7301** | **True** |

**Interpretation**: Ford-Roman inequality is violated at τ ≥ 1.0 for both massless and massive cases. This is expected for strong NEC violation (|ρ| ~ 15.6) and confirms the QFT caveat: the NEC-violating configuration is not protected by quantum inequalities at long sampling times.

### 3. Two-Way Coupling with Peak-Preserving Restrict (Exploratory)

**Test**: `v7_1_twoway_run.py` with `restrict_matter=True, matter_restrict_mode="peak", matter_blend=0.7`

| Step | ρ_min | n_neg | Phi_max | psi_max | Status |
|------|-------|-------|---------|---------|--------|
| 0 | -1.564e+01 | 27 | 3.000 | 2.828 | NEC violation |
| 25 | 0.0000 | 0 | 1.246 | 0.726 | **NEC lost** |
| 50 | 0.0000 | 0 | 0.254 | 0.000 | Phi decaying |
| 100 | 0.0000 | 0 | 0.132 | 0.000 | Phi → 0 |

**Conclusion**: Two-way coupling with peak-preserving restrict fails to maintain NEC violation. The restrict operation (even with peak preservation) causes Phi to decay to zero within ~25 steps. One-way coupling (coarse matter evolves independently, fine grid receives prolonged coarse matter) remains the canonical approach.

---

## Key Technical Decisions

### One-Way Coupling Confirmed Canonical
- **Why**: Two-way coupling loses NEC violation due to restrict averaging
- **How**: Coarse grid matter evolves independently; fine grid receives prolonged coarse matter; psi re-derived from `psi_equilibrium_smooth(Phi)` after prolongate
- **Result**: Stable NEC violation for 100+ steps on refined grid

### psi_equilibrium_smooth() with Smoothstep
- **Function**: `psi_equilibrium_smooth(Phi, m2, alpha, Eth, lam_psi, width=0.3)`
- **Purpose**: Avoids Laplacian spikes at Φ = Eth threshold
- **Used in**: v7.0 (complex+PTUM), v7.1 (Kerr+complex+PTUM)

### Kerr-Schild Metric with a=0.9
- **Source**: `kerr_schild_state()` from `bssn_ic_kerr_full.py`
- **Regularity**: |a| ≤ M satisfied (a=0.9, M=1.0)
- **Center**: α(0) ≈ 1, regular

---

## Unit Tests (All Passing)

| Test File | Tests | Status |
|-----------|-------|--------|
| `test_bssn_refined_kerr_complex_ptum.py` | 3 | ✅ |
| `test_bssn_refined_complex_ptum_v2.py` | 3 | ✅ |
| `test_bssn_refined_complex_ptum_v3.py` | 3 | ✅ |
| **Total** | **9** | **✅** |

---

## Files Modified

- `docs/roadmap.md` — v7.1 items marked ✅ complete
- `RELEASE_NOTES_v7.1-beta.md` — This file

---

## Git Actions

```bash
git add docs/roadmap.md RELEASE_NOTES_v7.1-beta.md
git commit -m "v7.1-beta consolidation: Kerr+complex+PTUM NEC stable, Ford-Roman violation documented, two-way coupling fails"
git tag v7.1-beta
git push origin main --tags
```

---

## Next: v7.2-beta

- Transient discharge attempt (negative result documented)
- Consolidate as v7.2.0 release