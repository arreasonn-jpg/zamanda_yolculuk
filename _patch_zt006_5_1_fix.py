"""Fix ZT-006.5.1: use properly normalized rest 4-velocity."""
from pathlib import Path

SRC = Path("src/zt005/run_zt006_5_1.py")
src = SRC.read_text(encoding="utf-8")

old = '''                sol = integrate_geodesic(
                    g, rv, zv,
                    x0=(0, 0, 0, 0),
                    U0=(1.0, 0, 0, 0),
                    tau_max=1.0,      # 1 second
                    n_steps=100,
                )'''

new = '''                # Properly normalized rest-frame 4-velocity:
                # g_tt U^t U^t = -1  =>  U^t = 1/sqrt(-g_tt)
                g_tt_origin = float(g[0, len(zv)//2, 0, 0])
                U_t0 = 1.0 / np.sqrt(max(-g_tt_origin, 1e-30))
                sol = integrate_geodesic(
                    g, rv, zv,
                    x0=(0, 0, 0, 0),
                    U0=(U_t0, 0, 0, 0),
                    tau_max=1.0,      # 1 second
                    n_steps=100,
                )'''

assert old in src
src = src.replace(old, new)
SRC.write_text(src, encoding="utf-8")
print("[OK] Fixed initial 4-velocity normalization")
print()
print("Now run:")
print("  python -m zt005.run_zt006_5_1 --mode quick")