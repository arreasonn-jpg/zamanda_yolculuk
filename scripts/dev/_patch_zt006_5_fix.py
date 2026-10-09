"""Fix ZT-006.5 — x0 must be 4D (t,x,y,z)."""
from pathlib import Path

TST = Path("tests/test_zt006_5.py")
SRC = Path("src/zt005/run_zt006_5.py")

# --- Fix source default ---
src = SRC.read_text(encoding="utf-8")
src = src.replace(
    "x0=(0.0, 0.0, 0.0),",
    "x0=(0.0, 0.0, 0.0, 0.0),"
)
SRC.write_text(src, encoding="utf-8")
print("[OK] Fixed src default x0")

# --- Fix tests ---
tst = TST.read_text(encoding="utf-8")
tst = tst.replace("x0=(0, 0, 0),", "x0=(0, 0, 0, 0),")
TST.write_text(tst, encoding="utf-8")
print("[OK] Fixed tests x0")

print()
print("Now run:")
print("  python -m pytest tests\\test_zt006_5.py -v")
print("  python -m zt005.run_zt006_5 --mode quick")