# v6.4-beta Release Notes

## Summary
**Full stack: Kerr + SC regular core + condensate + PTUM fields** — v6.4.0 release.

This consolidation completes the full integration of all three physics modules:
- **Kerr-Schild metric** with smoothed M_eff(r) regular core
- **SC condensate** (single real scalar, BCS-inspired potential)
- **PTUM fields** (fuel Φ + proxy ψ, tachyon condensation with NEC violation)

All configurations validated: single-grid (30-step, 200-step), refined 2:1 grid (100 steps), and strong NEC violation (A_fuel=3.0, 100 steps).

## Components Validated

### 1. Kerr-SC-PTUM Initial Conditions (`kerr_sc_ptum_ic.py`)
- Combines Kerr-SC metric + SC condensate + PTUM fields
- Returns full state dict with all 7 fields + potentials
- `combined_diagnostics()` computes total ρ, split ρ_SC/ρ_PTUM, negative cells

### 2. Single-Grid Driver (`bssn_rk4_sc_ptum.py`)
- RK4 evolution coupling BSSN + SC condensate + PTUM
- Returns 7-tuple: (state, φ_sc, Π_sc, Φ, Π_Φ, ψ, Π_ψ)
- Handles all stress-energy sources

### 3. Refined 2:1 Grid Driver (`bssn_refined_sc_ptum.py`)
- Coarse (N_c=9, L_c=20) + Fine (N_f=9, L_f=10) grids
- Prolongation: coarse → fine via tensor refinement
- Restriction: fine → coarse via averaging
- Full cycle: coarse RK4 → prolongate → fine RK4 → restrict

## Test Results

### Unit Tests (26 tests passing)
- `test_kerr_sc_ic.py`: 9 tests
- `test_sc_condensate.py`: 5 tests
- `test_kerr_sc_ptum_ic.py`: 6 tests
- `test_bssn_refined_sc_ptum.py`: 3 tests
- `test_bssn_rk4_sc_ptum.py`: 3 tests

### Single-Grid Evolution (30 steps, a=0.0/0.5/0.9)
```
FULL SC+PTUM EVOLUTION: N=9 dt=0.0001 steps=30
a=0.00 rho_min=-2.842e-01 rho_sc=6.553e-150 rho_pt=-2.842e-01 Kmax=3.628e-01 alpha_min=0.8776 fin=True
a=0.50 rho_min=-2.842e-01 rho_sc=6.553e-150 rho_pt=-2.842e-01 Kmax=3.584e-01 alpha_min=0.8766 fin=True
a=0.90 rho_min=-2.842e-01 rho_sc=6.553e-150 rho_pt=-2.842e-01 Kmax=3.497e-01 alpha_min=0.8742 fin=True
```
- NEC violation maintained (ρ_min ≈ -0.284)
- SC contribution negligible (ρ_sc ~ 10⁻¹⁵⁰)
- All spins stable, finite fields

### Single-Grid Long Run (200 steps, a=0.9)
```
FULL STACK LONG: N=9 a=0.9 dt=0.0001 steps=200
step=  40 rho_min=-2.8422e-01 Kmax=3.5695e-01 alpha_min=0.8742 fin=True
step=  80 rho_min=-2.8423e-01 Kmax=3.8583e-01 alpha_min=0.8742 fin=True
step= 120 rho_min=-2.8424e-01 Kmax=4.1477e-01 alpha_min=0.8742 fin=True
step= 160 rho_min=-2.8425e-01 Kmax=4.4378e-01 alpha_min=0.8742 fin=True
step= 200 rho_min=-2.8427e-01 Kmax=4.7287e-01 alpha_min=0.8742 fin=True
```
- NEC violation stable at ρ ≈ -0.284
- K_max grows from 0.36 → 0.47 (controlled)
- α_min constant at 0.8742 (no collapse)
- All fields finite

### Refined 2:1 Grid (100 steps, A_fuel=1.5)
```
REFINED FULL STACK: N_c=9 N_f=9 dt=0.0002 steps=100
step=  25 rho_min=8.0934e-10 Kmax_c=8.043e-02 Kmax_f=8.032e-02 alpha_min=0.8805 fin=True
step=  50 rho_min=8.0934e-10 Kmax_c=7.783e-02 Kmax_f=7.773e-02 alpha_min=0.8805 fin=True
step=  75 rho_min=8.0934e-10 Kmax_c=7.524e-02 Kmax_f=7.514e-02 alpha_min=0.8805 fin=True
step= 100 rho_min=8.0934e-10 Kmax_c=7.265e-02 Kmax_f=7.255e-02 alpha_min=0.8805 fin=True
```
- K_max decreasing on both grids (7.8% → 7.3%)
- ρ_min ~ 0 (no NEC violation at A_fuel=1.5)
- Coarse/fine K_max track closely

### Refined 2:1 Grid Strong (100 steps, A_fuel=3.0, NEC violation)
```
REFINED FULL STACK (STRONG): A_fuel=3.0 N_c=9 N_f=9
initial rho_total_min=-1.564e+01 rho_ptum_min=-1.564e+01
step=  25 rho_min=-9.8550e-02 rho_pt=-9.8550e-02 n_neg=3 Kmax_c=9.486e-02 Kmax_f=9.534e-02 alpha_min=0.8805 fin=True
step=  50 rho_min=-9.8553e-02 rho_pt=-9.8553e-02 n_neg=3 Kmax_c=1.067e-01 Kmax_f=1.072e-01 alpha_min=0.8805 fin=True
step=  75 rho_min=-9.8559e-02 rho_pt=-9.8559e-02 n_neg=3 Kmax_c=1.186e-01 Kmax_f=1.191e-01 alpha_min=0.8805 fin=True
step= 100 rho_min=-9.8566e-02 rho_pt=-9.8566e-02 n_neg=3 Kmax_c=1.305e-01 Kmax_f=1.310e-01 alpha_min=0.8805 fin=True
```
- NEC violation maintained (ρ_min ≈ -0.099, 3 negative cells)
- K_max grows from 0.095 → 0.131 (controlled)
- Coarse/fine track closely
- All fields finite

## Files in This Release

### Core Implementation
- `src/zt005/kerr_sc_ptum_ic.py` — Full stack IC (metric + SC + PTUM)
- `src/zt005/bssn_rk4_sc_ptum.py` — Single-grid RK4 driver
- `src/zt005/bssn_refined_sc_ptum.py` — Refined 2:1 grid driver
- `src/zt005/matter/sc_condensate.py` — SC condensate
- `src/zt005/matter/ptum_field.py` — PTUM fuel/proxy fields
- `src/zt005/matter/matter_scalar.py` — Scalar stress-energy
- `src/zt005/matter/matter_ptum.py` — PTUM stress-energy

### Tests
- `tests/test_kerr_sc_ic.py` — 9 tests
- `tests/test_sc_condensate.py` — 5 tests
- `tests/test_kerr_sc_ptum_ic.py` — 6 tests
- `tests/test_bssn_refined_sc_ptum.py` — 3 tests
- `tests/test_bssn_rk4_sc_ptum.py` — 3 tests

### Validation Scripts
- `scripts/v6_sc_ptum_full_run.py` — 30-step parameter sweep
- `scripts/v6_sc_ptum_full_long.py` — 200-step long run
- `scripts/v6_refined_sc_ptum_run.py` — 100-step refined (A_fuel=1.5)
- `scripts/v6_refined_sc_ptum_run_strong.py` — 100-step refined (A_fuel=3.0)

## Next Steps (v6.5)
- Refined full stack 200+ step long run
- Complex scalar + PTUM with equilibrium seed
- Ford-Roman quantum inequality on full stack

## Git Tag
```bash
git tag v6.4-beta
git push origin v6.4-beta
```