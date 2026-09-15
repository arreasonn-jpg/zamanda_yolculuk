# ZT-006.2.17
Metric convention / weak-field sanity gate.

ZT-006.2.16 produced impossible determinant/timelike diagnostics for a
perturbation of order 1e-46. This stage validates the perturbation directly
with g = eta + sym(h)/2 and compares raw Green output against the trace-reversed
form.

Run:
python -m pytest
python -m zt005.run_zt006_2_17 --mode quick
python -m zt005.run_zt006_2_17 --mode full
