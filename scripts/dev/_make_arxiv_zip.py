"""Create arXiv-compatible ZIP (forward slashes)."""
import zipfile
from pathlib import Path

PAPER_DIR = Path("paper")
OUTPUT_ZIP = Path("arxiv_submission_v2.zip")

# Files to include (relative to paper/)
FILES = [
    "main.tex",
    "refs.bib",
    "figures/fig1_unified_mechanisms.png",
    "figures/fig2_em_scaling.png",
    "figures/fig3_best_case_tiers.png",
    "figures/fig4_monte_carlo.png",
    "figures/fig5_sensitivity_tornado.png",
    "figures/fig6_backpack_layout.png",
]

if OUTPUT_ZIP.exists():
    OUTPUT_ZIP.unlink()

with zipfile.ZipFile(OUTPUT_ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
    for rel in FILES:
        src = PAPER_DIR / rel
        if not src.exists():
            print(f"[MISSING] {src}")
            continue
        # Force forward-slash arcname
        arcname = rel.replace("\\", "/")
        zf.write(src, arcname=arcname)
        print(f"[OK] {arcname}")

# Verify structure
print()
print("=== ZIP contents ===")
with zipfile.ZipFile(OUTPUT_ZIP, "r") as zf:
    for info in zf.namelist():
        print(f"  {info}")
        assert "\\" not in info, f"BACKSLASH in {info}!"
print()
print(f"[OK] {OUTPUT_ZIP} created, size = {OUTPUT_ZIP.stat().st_size / 1024:.1f} KB")
print("[OK] All paths use forward slashes")