"""Fix ZT-007B.4 EM gap test."""
from pathlib import Path

t = Path("tests/test_zt007_b4.py")
c = t.read_text(encoding="utf-8")

old = '''def test_em_gap_exceeds_1e40():
    d = collect_mechanism_data()
    assert d["EM"]["gap"] > 1e40'''

new = '''def test_em_gap_exceeds_expected():
    """EM current gap is ~10^20; energy-density gap is ~10^40."""
    d = collect_mechanism_data()
    assert d["EM"]["gap"] > 1e19
    assert d["EM"]["rho_gap"] > 1e40'''

assert old in c, "old test block not found"
c = c.replace(old, new)
t.write_text(c, encoding="utf-8")
print("[OK] Fixed test_zt007_b4.py")