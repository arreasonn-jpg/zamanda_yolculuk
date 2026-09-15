# ZT-006.2 — Source Convergence + Current Scaling

This stage validates the weak-field source calculation before adding Earth/world sources.

Source convergence: 200, 400, 800, 1600 nested deterministic Monte-Carlo source samples.
Current scaling: 25, 50, 100, 200, 400 A.

Expected in the linear EM + linearized-GR regime:
B ∝ I, E ∝ I, T ∝ I², h ∝ I².

Run:
```powershell
pip install -e .
pytest
python -m zt005.run_zt006_2
```
For larger/faster runs use WSL/Linux.
No CTC conclusion is made here.
