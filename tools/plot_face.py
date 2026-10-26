"""Plot the top-face height map written by `urea-kmc face` (like the paper's AFM image, Fig. 1c).

Run: uv run python tools/plot_face.py runs/face-methanol     -> runs/face-methanol/face.png
"""
import sys, pathlib
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "runs/face-water")
h = np.loadtxt(run / "top.csv", delimiter=",", ndmin=2)
along = h.shape[1] * 2 * 5.565 / np.sqrt(2) / 10  # nm along [1-10]
up = h.shape[0] * 0.4684                           # nm along c
log = [l for l in (run / "log.csv").read_text().splitlines() if l and l[0].isdigit()]
grown = log[-1].split(",")[3] if log else "?"
fig, ax = plt.subplots(figsize=(7, 7 * up / along + 0.8))
vmax = max(1, int(np.nanmax(h)))
im = ax.imshow(h, origin="lower", extent=[0, along, 0, up], cmap="viridis", vmin=min(0, int(np.nanmin(h))), vmax=vmax,
               interpolation="nearest")
ax.set_xlabel("nm along [1-10]")
ax.set_ylabel("nm along c")
ax.set_title(f"{run.name}: [110] face, height in layers\n(mean growth {grown} layers per face)", fontsize=10)
fig.colorbar(im, ax=ax, label="layers (0.39 nm each)", shrink=0.8)
fig.tight_layout()
out = run / "face.png"
fig.savefig(out, dpi=130)
print(f"wrote {out}")
