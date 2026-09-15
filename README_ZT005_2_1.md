# ZT-005.2.1 — Optimized Wire-Distance Diagnostics

This stage separates grid-resolution study from filament resolution and
diagnoses how close active-volume sample points are to the idealized
zero-thickness source wires.

The distance computation uses `scipy.spatial.cKDTree` over a dense sampling
of each source polyline, avoiding the O(N_points * N_segments) brute-force
scan that made the previous version too slow.

Diagnostic thresholds:
- 2 dx
- 3 dx
- 4 dx

These are not physical exclusions. They quantify regions where a
zero-thickness filament representation may be numerically poorly conditioned.

Run:
    pip install -e .
    pytest
    python -m zt005.run_zt005_2_1
