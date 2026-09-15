# Scripts

Runnable entry points currently live in the package itself as
`src/zt005/run_*.py` modules, invoked as:

```bash
python -m zt005.run_zt006_3_gr_layer --mode quick   # GR validation layer
python -m zt005.run_zt006_2_20 --mode quick          # historical stage
python -m zt005.run_zt005                            # ZT-005 EM consistency
```

Output policy:

* Historical stage runners (re-runs) write to `results/exploratory/`.
  `results/archived/` holds the committed historical record and is never
  overwritten by runners.
* The ZT-006.3 GR layer writes to `results/validated/`.
* Ad-hoc scans belong in `results/exploratory/` and are never cited.

This directory is reserved for future operational scripts (CI helpers,
result promotion/archival tooling).
