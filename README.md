# ZAMANDA YOLCULUK — ZT-005
## Maxwell-Consistent Magnetoquasistatic Field Validation

ZT-005 replaces the ZT-004 ad-hoc electric-field proxy with a vector-potential
based magnetoquasistatic model.

For a closed filamentary current distribution:

A(r) = mu0/(4*pi) ∫ I dl / |r-r'|

For harmonic current I(t)=I0 cos(omega t):

A(r,t)=A0(r) cos(omega t)

B(r,t)=curl A0(r) cos(omega t)

E(r,t)=omega A0(r) sin(omega t)

under the magnetoquasistatic Coulomb-gauge baseline for closed current loops.

The code then performs numerical consistency checks:

1. div(B) ≈ 0
2. curl(E) ≈ -dB/dt
3. B from curl(A) versus Biot-Savart B

### Scientific scope

This is still not a full-wave finite-element Maxwell solver.
It is a controlled magnetoquasistatic validation stage.

The purpose is to determine whether the electromagnetic field model used in
ZT-004 is internally consistent before using its stress-energy tensor in
relativistic calculations.

No CTC, time machine, or macroscopic spacetime manipulation is assumed.

### Run

```powershell
pip install -e .
pytest
python -m zt005.run_zt005
```

Outputs:
- `zt005_results.json`
- `zt005_results/*.png`
