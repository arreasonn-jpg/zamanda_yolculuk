"""scripts/dev/refactor_v13.py — one-shot subpackage split (v13-alpha).

Moves 14 core modules into two subpackages and rewrites imports
repo-wide in one pass:

  amr/    : amr_interp, amr_bh_track, refinement, refinement_tensor,
            bssn_amr_consistent
  matter/ : matter_scalar, matter_complex, matter_ptum, matter_proca,
            matter_brans_dicke, sc_condensate, exotic_matter,
            phase_transition_common, ptum_field

Rules:
- absolute  zt005.<M>(...)         -> zt005.<pkg>.<M>
- relative  from .<M> import       -> from .<pkg>.<M> import
  (both applied ONLY to files outside the moved subpackages, so
   intra-subpackage relative imports stay valid)
- files are moved with `git mv` (history preserved)
- dry run prints the plan; pass --apply to execute
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(".")
GROUPS = {
    "amr": ["amr_interp", "amr_bh_track", "refinement", "refinement_tensor",
            "bssn_amr_consistent"],
    "matter": ["matter_scalar", "matter_complex", "matter_ptum",
               "matter_proca", "matter_brans_dicke", "sc_condensate",
               "exotic_matter", "phase_transition_common", "ptum_field"],
}
SCAN_DIRS = ["src", "tests", "scripts", "stages"]  # md: none reference code
# scan only python files; stages has no python (kept for future)
APPLY = "--apply" in sys.argv


def iter_target_files():
    for d in SCAN_DIRS:
        dpath = ROOT / d
        if not dpath.exists():
            continue
        for p in dpath.rglob("*.py"):
            s = str(p)
            if "__pycache__" in s or "egg-info" in s:
                continue
            if "/src/zt005/amr/" in s or "/src/zt005/matter/" in s:
                continue
            yield p


def main():
    # 1) import rewrites (before moving: rewrite then git mv)
    changed = {}
    for p in iter_target_files():
        text = p.read_text(errors="ignore")
        orig = text
        for pkg, mods in GROUPS.items():
            for m in mods:
                text = text.replace(f"zt005.{m}", f"zt005.{pkg}.{m}")
                text = re.sub(rf"from \.{re.escape(m)} import",
                              f"from .{pkg}.{m} import", text)
                text = re.sub(rf"from \.{re.escape(m)} \(", 
                              f"from .{pkg}.{m} (", text)
        if text != orig:
            changed[p] = (orig.count("\n"), text.count("\n"))
            if APPLY:
                p.write_text(text)
    print(f"files with imports rewritten: {len(changed)}")

    # 2) moves
    for pkg, mods in GROUPS.items():
        pkg_dir = ROOT / "src" / "zt005" / pkg
        if APPLY:
            pkg_dir.mkdir(exist_ok=True)
            (pkg_dir / "__init__.py").write_text("")
        for m in mods:
            src = ROOT / "src" / "zt005" / f"{m}.py"
            dst = pkg_dir / f"{m}.py"
            print(("MOVE " if APPLY else "WOULD MOVE "),
                  src, "->", dst, "(exists)" if dst.exists() else "")
            if APPLY and src.exists():
                subprocess.run(["git", "mv", str(src), str(dst)], check=True)


if __name__ == "__main__":
    main()
