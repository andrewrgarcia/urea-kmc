"""Plot a urea-kmc surface dump: three projections, shaded by coordination (Fig. 4 of the paper:
dark, medium, light blue = 3, 4, 5 neighbors).

Run: uv run python tools/plot.py runs/water        -> runs/water/morphology.png
"""
import sys, pathlib
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "runs/water")
lines = (run / "surface.csv").read_text().splitlines()
if len(lines) < 2:
    print(f"{run}/surface.csv is empty: the crystal dissolved, nothing to plot")
    sys.exit(0)
d = np.loadtxt(run / "surface.csv", delimiter=",", skiprows=1, ndmin=2)
x, y, z, n = d[:, 0], d[:, 1], d[:, 2], d[:, 4].astype(int)
u, v = (x + y) / np.sqrt(2), (x - y) / np.sqrt(2)  # along [110] and [1-10]
shade = {0: "#dddddd", 1: "#dddddd", 2: "#bbbbbb", 3: "#08306b", 4: "#2171b5", 5: "#9ecae1"}
colors = [shade[k] for k in n]

fig, ax = plt.subplots(1, 3, figsize=(12, 6))
for a, (h, vv, title) in zip(ax, [(u, z, "side, along [1-10]"), (v, z, "side, along [110]"), (u, v, "top, along c")]):
    order = np.argsort(n)[::-1]  # draw high coordination first, kinks on top
    a.scatter(np.asarray(h)[order], np.asarray(vv)[order], s=4, c=np.asarray(colors)[order], linewidths=0)
    a.set_aspect("equal")
    a.set_title(title)
    a.set_xlabel("nm")
ax[0].set_ylabel("nm")
fig.suptitle(f"{run.name}: {len(n)} surface molecules")
fig.tight_layout()
out = run / "morphology.png"
fig.savefig(out, dpi=150)
print(f"wrote {out}")
