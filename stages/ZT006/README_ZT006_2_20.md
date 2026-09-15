# ZT-006.2.20 — External-Observer Green-Kernel Convergence

ZT-006.2.19 correctly failed: the final-grid Green result was not converged
for observers inside the source volume.

This stage is an independent numerical control. Observation points are all
outside the active cylinder, so the Green kernel 1/|x-x'| has no singularity.

Run:
```powershell
python -m pip install -e .
python -m pytest
python -m zt005.run_zt006_2_20 --mode quick
python -m zt005.run_zt006_2_20 --mode full
```

No CTC conclusion is made here.
