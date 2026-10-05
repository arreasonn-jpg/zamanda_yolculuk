"""Fix ZT-006.5 main() x0 call."""
from pathlib import Path

SRC = Path("src/zt005/run_zt006_5.py")
src = SRC.read_text(encoding="utf-8")

# Replace all remaining 3D x0 calls with 4D
n1 = src.count("x0=(0, 0, 0),")
src = src.replace("x0=(0, 0, 0),", "x0=(0, 0, 0, 0),")

n2 = src.count("x0=(0, 0, 0)")
if n2 > 0:
    src = src.replace("x0=(0, 0, 0)", "x0=(0, 0, 0, 0)")

SRC.write_text(src, encoding="utf-8")
print(f"[OK] Replaced {n1 + n2} x0 call(s)")
print()
print("Now run:")
print("  python -m pytest tests\\test_zt006_5.py -v")
print("  python -m zt005.run_zt006_5 --mode quick")