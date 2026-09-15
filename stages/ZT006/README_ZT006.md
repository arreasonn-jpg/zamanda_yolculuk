
# ZT-006 — Linearized Einstein / Gravitational-Source Diagnostics

This stage takes the validated pointwise electromagnetic stress-energy and
computes controlled weak-field source diagnostics:

- T00, T0i, Tij
- trace-reversed source
- characteristic h00 and h0i scales
- proper-time fractional-shift proxy

It does NOT solve the full Einstein equations and does NOT claim a CTC.

Earth and device background stress-energy are intentionally omitted until a
separate explicit world model is defined.

Run:
pip install -e .
pytest
python -m zt005.run_zt006
