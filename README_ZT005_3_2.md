# ZT-005.3.2 — Sparse / Pointwise Cross-Validation

Uses a finite set of points inside the active cylinder instead of constructing
a full 3-D vector-potential tensor.

For each geometry:
- Biot-Savart B is evaluated at the validation points.
- curl(A) is evaluated pointwise with 2nd- and 4th-order stencils.
- Points closer than 3h to a source wire are omitted for that h.

Step sizes:
- h = 0.02 m
- h = 0.01 m
- h = 0.005 m

Run:
```powershell
pip install -e .
pytest
python -m zt005.run_zt005_3_2
```

No GR/CTC conclusion is made here.
