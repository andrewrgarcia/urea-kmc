"""Paper figures from the runs written by `makers figure-runs`.

Run: uv run python tools/figures.py [RUNS_DIR] [OUT_DIR]
     defaults: runs/fig  ../../papers/urea-kmc-replication/figures  (both relative to this repo)

Writes morphology.pdf, faces.pdf, island.pdf.
"""
import sys, pathlib
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle

root = pathlib.Path(__file__).resolve().parent.parent
runs = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else root / "runs/fig"
out = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else root.parent.parent / "papers/urea-kmc-replication/figures"
out.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "figure.dpi": 150})
SHADE = {3: "#08306b", 4: "#2171b5", 5: "#9ecae1"}  # paper's Fig. 4: dark, medium, light blue = 3, 4, 5 neighbors
A, C = 0.5565, 0.4684  # nm
D110 = A / np.sqrt(2)


def surface(run):
    d = np.loadtxt(runs / run / "surface.csv", delimiter=",", skiprows=1, ndmin=2)
    x, y, z, n = d[:, 0], d[:, 1], d[:, 2], d[:, 4].astype(int)
    u, v = (x + y) / np.sqrt(2), (x - y) / np.sqrt(2)
    keep = (n >= 3) & (n <= 5)
    return u[keep] - u[keep].mean(), v[keep] - v[keep].mean(), z[keep] - z[keep].min(), n[keep]


def scatter(ax, h, vv, n):
    order = np.argsort(n)[::-1]  # five-coordinated first, kinks on top
    ax.scatter(h[order], vv[order], s=0.6, c=[SHADE[k] for k in n[order]], linewidths=0, rasterized=True)
    ax.set_aspect("equal")


# ---- Figure: morphologies (two side views each, so the 90-degree twist of the ends shows) ------------------------
fig, axes = plt.subplots(1, 4, figsize=(7.0, 4.6), sharey=True)
panels = [("water", "water, C = 0.24", 0), ("water", "water, C = 0.24", 1),
          ("methanol", "methanol, C = 0.21", 0), ("methanol", "methanol, C = 0.21", 1)]
cache = {}
for ax, (run, label, view) in zip(axes, panels):
    if run not in cache:
        cache[run] = surface(run)
    u, v, z, n = cache[run]
    scatter(ax, u if view == 0 else v, z, n)
    ax.set_title(label + "\n" + ("viewed along [1$\\bar{1}$0]" if view == 0 else "viewed along [110]"))
    ax.set_xlabel("nm")
    ax.set_xticks([-10, 0, 10])
axes[0].set_ylabel("nm along c")
handles = [plt.Line2D([], [], marker="o", ls="", color=SHADE[k], label=f"{k} neighbors") for k in (3, 4, 5)]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig(out / "morphology.pdf", dpi=200)
plt.close(fig)


# ---- Figure: {110} face height maps -----------------------------------------------------------------------------
def heights(run):
    return np.loadtxt(runs / run / "top.csv", delimiter=",", ndmin=2)


def show_map(ax, h, title):
    along, up = h.shape[1] * 2 * D110, h.shape[0] * C
    lo, hi = int(min(0, h.min())), int(max(1, h.max()))
    levels = np.arange(lo, hi + 1)
    cmap = ListedColormap(plt.cm.viridis(np.linspace(0, 1, len(levels))))
    norm = BoundaryNorm(np.arange(lo - 0.5, hi + 1.5), cmap.N)
    im = ax.imshow(h, origin="lower", extent=[0, along, 0, up], cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("nm along [1$\\bar{1}$0]")
    return im, levels


fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.9))
for ax, (run, title) in zip(axes, [("face-water-0.24", "water, C = 0.24"),
                                   ("face-methanol-0.21", "methanol, C = 0.21"),
                                   ("face-methanol-0.19", "methanol, C = 0.19")]):
    im, levels = show_map(ax, heights(run), title)
    cb = fig.colorbar(im, ax=ax, shrink=0.75, ticks=levels)
    cb.set_label("layers")
axes[0].set_ylabel("nm along c")
fig.tight_layout()
fig.savefig(out / "faces.pdf", dpi=200)
plt.close(fig)


# ---- Figure: seeded island ----------------------------------------------------------------------------------------
isl = runs / "island-methanol-0.19"
h = heights("island-methanol-0.19")
d = np.loadtxt(isl / "island.csv", delimiter=",", skiprows=2, ndmin=2)
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1, 1.3]})
im, levels = show_map(ax0, h, "methanol, C = 0.19, 100 ns, top face")
ax0.set_ylabel("nm along c")
period, nz = h.shape[1], h.shape[0]
wi, hi = 25, 43  # the seeded island: 25 columns x 43 rows (20 x 20 nm)
ax0.add_patch(Rectangle(((period - wi) // 2 * 2 * D110, (nz - hi) // 2 * C), wi * 2 * D110, hi * C,
                        fill=False, ec="white", lw=0.8, ls="--"))
fig.colorbar(im, ax=ax0, shrink=0.75, ticks=levels).set_label("layers")
t = d[:, 0]
for col, face, ls in [(1, "top", "-"), (5, "bottom", "--")]:
    ax1.plot(t, d[:, col], ls, color="#c2410c", label=f"{face} face, toward +c")
    ax1.plot(t, d[:, col + 1], ls, color="#1d4ed8", label=f"{face} face, toward \N{MINUS SIGN}c")
ax1.set_xlabel("time (ns)")
ax1.set_ylabel("edge displacement (nm)")
ax1.legend(frameon=False)
ax1.set_title("island edges")
fig.tight_layout()
fig.savefig(out / "island.pdf", dpi=200)
plt.close(fig)

print(f"wrote {out}/morphology.pdf, faces.pdf, island.pdf")
