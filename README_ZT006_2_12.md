# ZT-006.2.12 — Current-Conserving Finite-Cross-Section Model

ZT-006.2.11 accidentally applied the full parent current to every subfilament.
That artificially multiplied the source strength by the number of
cross-section filaments.

This version enforces:

    I_each = I_total / N_filaments

and compares finite-radius results against the zero-radius baseline at fixed
validation points.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_12 --mode quick
python -m zt005.run_zt006_2_12 --mode full
```
