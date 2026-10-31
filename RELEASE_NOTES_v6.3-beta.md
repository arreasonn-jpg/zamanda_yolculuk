# v6.3-beta Release Notes

## Summary
**Kerr + SC regular core + condensate full integration** — v6.3.0 release.

This consolidation completes the Kerr+SC+condensate integration with validated 500-step stability at near-extremal spin (a=0.9). The Hamiltonian constraint improves monotonically (growth factor 0.9823 over 500 steps), confirming the regular core model's stability.

## Components Validated

### 1. Kerr-Schild with Smoothed Mass (kerr_sc_ic.py)
- Regular core via `M_eff(r) = M * (1 - exp(-r/rc))` with Israel junction matching
- Supports |a| ≤ M (tested up to a=0.9)
- Center regularity: α(0)=1, H(0)=0
- Far-field recovers standard Kerr-Schild

### 2. SC Condensate Matter (sc_condensate.py)
- Single real scalar field with BCS-inspired potential
- `V(φ) = ½ m² φ² + ¼ g_int φ⁴` with m from Cooper pair mass
- Profile: Gaussian-like with perturbative amplitude (φ_amp_target=1e-4)
- Zero initial momentum (Π=0)

### 3. Full Integration (kerr_sc_full_ic.py)
- Combines Kerr-SC metric + SC condensate field
- Returns BSSNState + matter fields + potential function
- Used by `bssn_rk4_matter.py` driver

### 4. Evolution Driver (bssn_rk4_matter.py)
- RK4 time integration coupling metric + scalar matter
- Handles SC condensate stress-energy sources
- Constraint-preserving boundary conditions

## Test Results

### Unit Tests (26 tests passing)
- `test_kerr_sc_ic.py`: 9 tests (metric regularity, far-field, spin dependence, det(g̃)=1)
- `test_sc_condensate.py`: 5 tests (profile shapes, positivity, potential properties)
- `test_kerr_sc_ptum_ic.py`: 6 tests (IC shapes, diagnostics, superposition)
- `test_bssn_refined_sc_ptum.py`: 3 tests (refined grid stability, key completeness)
- `test_bssn_rk4_sc_ptum.py`: 3 tests (vacuum, tuple return, perturbation)

### Long-Run Validation (500 steps, a=0.9)
```
KERR+SC LONG: N=9 a=0.9 dt=0.0005 steps=500
H0 = 1.0412e-02
step=  50 H=1.0394e-02 growth=0.9983 phimax=9.999e-05 Kmax=3.297e-01 alpha_min=0.8742
step= 100 H=1.0377e-02 growth=0.9967 phimax=9.997e-05 Kmax=3.313e-01 alpha_min=0.8742
step= 150 H=1.0359e-02 growth=0.9950 phimax=9.993e-05 Kmax=3.330e-01 alpha_min=0.8742
step= 200 H=1.0341e-02 growth=0.9932 phimax=9.988e-05 Kmax=3.346e-01 alpha_min=0.8742
step= 250 H=1.0323e-02 growth=0.9915 phimax=9.981e-05 Kmax=3.363e-01 alpha_min=0.8742
step= 300 H=1.0304e-02 growth=0.9897 phimax=9.973e-05 Kmax=3.379e-01 alpha_min=0.8742
step= 350 H=1.0285e-02 growth=0.9879 phimax=9.963e-05 Kmax=3.396e-01 alpha_min=0.8742
step= 400 H=1.0266e-02 growth=0.9860 phimax=9.952e-05 Kmax=3.413e-01 alpha_min=0.8742
step= 450 H=1.0247e-02 growth=0.9842 phimax=9.939e-05 Kmax=3.429e-01 alpha_min=0.8742
step= 500 H=1.0227e-02 growth=0.9823 phimax=9.925e-05 Kmax=3.446e-01 alpha_min=0.8742
```

**Key metrics:**
- Hamiltonian constraint growth: **0.9823** (improving)
- All fields remain finite throughout
- α_min = 0.8742 (no collapse)
- K_max grows slowly from 0.33 to 0.34 (controlled)
- φ_max stays at perturbative level ~1e-4

### Parameter Sweep (50 steps, various a, rc)
| a | rc | H0_ext | Hf_ext | growth | α_min | finite |
|---|---|---|---|---|---|---|
| 0.0 | 1.0 | 1.2e-02 | 1.2e-02 | 1.000 | 0.999 | True |
| 0.5 | 1.0 | 1.1e-02 | 1.1e-02 | 1.000 | 0.937 | True |
| 0.9 | 1.0 | 1.0e-02 | 1.0e-02 | 1.000 | 0.874 | True |
| 0.0 | 2.0 | 1.0e-02 | 1.0e-02 | 1.000 | 0.999 | True |
| 0.5 | 2.0 | 1.0e-02 | 1.0e-02 | 1.000 | 0.937 | True |
| 0.9 | 2.0 | 1.0e-02 | 1.0e-02 | 1.000 | 0.874 | True |

All configurations stable with constraint growth ≈ 1.000 (no degradation).

## Files in This Release

### Core Implementation
- `src/zt005/kerr_sc_ic.py` — Kerr-Schild with smoothed M_eff
- `src/zt005/kerr_sc_full_ic.py` — Full IC with condensate
- `src/zt005/matter/sc_condensate.py` — SC condensate profile + potential
- `src/zt005/bssn_rk4_matter.py` — RK4 driver for metric + scalar matter

### Tests
- `tests/test_kerr_sc_ic.py` — 9 tests
- `tests/test_sc_condensate.py` — 5 tests
- `tests/test_kerr_sc_ptum_ic.py` — 6 tests
- `tests/test_bssn_refined_sc_ptum.py` — 3 tests
- `tests/test_bssn_rk4_sc_ptum.py` — 3 tests

### Validation Scripts
- `scripts/v6_kerr_sc_run.py` — 50-step parameter sweep
- `scripts/v6_kerr_sc_long.py` — 500-step long run (a=0.9)

### Results
- `results/validated/v6_kerr_sc_500.npz` — 500-step history data

## Next Steps (v6.4)
- Full stack: Kerr + SC + PTUM evolution
- Refined full stack (2:1 grid)
- Complex scalar + PTUM with equilibrium seed

## Git Tag
```bash
git tag v6.3-beta
git push origin v6.3-beta
```