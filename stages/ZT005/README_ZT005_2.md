
# ZT-005.2 — Production-Resolution Field Convergence

ZT-005.2 returns to G1–G4 after the analytic-loop benchmark of ZT-005.1.

It tests whether the large ZT-005 discrepancy between `curl(A)` and
Biot–Savart decreases as both spatial and filament resolution increase.

## Default study

```powershell
pip install -e .
pytest
python -m zt005.run_zt005_2
```

Default resolutions:

- grid side: 25, 41, 61
- y-grid: scaled to the 2.5 m cylinder height
- circular filament: 128, 256, 512 segments
- helical filament: 240, 480 segments (paired with the circular settings)

The study is intentionally diagnostic. A geometry is not declared "best"
from one resolution. The intended criterion is convergence of the field
metrics.

For a faster smoke run:

```powershell
python -m zt005.run_zt005_2 --grids 17,25 --circle-segments 64,128 --helix-segments 96,192
```

No GR, CTC, causality, or time-travel claim is made at this stage.
