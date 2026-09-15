# ZT-006.2.13 — Parent-Conductor Current-Conserving Finite-Wire Model

Important correction: G1/G2/G3/G4 contain multiple parent conductors. The
configured current I belongs to EACH parent conductor. When a parent conductor
is thickened, only its cross-section subfilaments divide that parent's I.

Thus:
    sum(I_subfilaments within one parent) = I_parent

and, at radius=0, each parent is exactly one filament carrying I_parent.

Run:
```powershell
pip install -e .
python -m pytest
python -m zt005.run_zt006_2_13 --mode quick
python -m zt005.run_zt006_2_13 --mode full
```
