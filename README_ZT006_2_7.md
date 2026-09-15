# ZT-006.2.7 — Stress-Energy / Field Component Conditioning

ZT-006.2.6 showed that G4 is not uniquely associated with unusually high
normalized field gradients. This stage therefore asks whether the large
full-tensor sensitivity previously seen for G4 is caused by cancellation or
conditioning inside T_mu_nu itself.

For h = 0.02, 0.01, 0.005, 0.0025, compare:
- E relative error
- B relative error
- T00 relative RMS
- T0i relative RMS
- Tij relative RMS

The 0.005 result is the reference.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_7 --mode quick
python -m zt005.run_zt006_2_7 --mode full
```
No GR/CTC conclusion is made here.
