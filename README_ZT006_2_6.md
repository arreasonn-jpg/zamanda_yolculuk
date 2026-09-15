# ZT-006.2.6 — Field Gradient / Local Conditioning

Tests whether remaining geometry instability is associated with high local field
gradients rather than wire proximity or source cancellation.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_6 --mode quick
python -m zt005.run_zt006_2_6 --mode full
```

No GR/CTC conclusion is made here.
