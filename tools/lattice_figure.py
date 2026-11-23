"""Schematic of the six neighbors of a urea molecule, drawn from the same coordinates the code uses.

Run: uv run python tools/lattice_figure.py [OUT_DIR]   (default: ../../papers/urea-kmc-replication/figures)
Writes lattice-check.pdf, the reference for the hand-drawn Fig. 1 (figures/urea-neighbors.svg -> lattice.pdf,
which this script must not overwrite): left, the molecule's own plane (its two type-1 partners, and A/B along c);
right, the view along c (all four in-plane partners around the central molecule).
Labels follow the paper's Fig. 2 as drawn (A on the NH2 side); the table rows are read the other way (Sec. 3.1).
"""
import sys, pathlib
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = pathlib.Path(__file__).resolve().parent.parent
out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else root.parent.parent / "papers/urea-kmc-replication/figures"
out.mkdir(parents=True, exist_ok=True)

A, C, ZC = 5.565, 4.684, 0.3260  # Angstrom, fractional z of carbon (lattice.rs)
# alpha molecule at (0, 1/2, zC), C=O toward +c; beta at (1/2, 0, 1-zC), C=O toward -c (lattice.rs)
def pos(s, i, j, k):
    return np.array([i, j + 0.5, k + ZC]) if s == 0 else np.array([i + 0.5, j, k + 1 - ZC])

# Central molecule: the beta molecule of cell (0,0,0), whose C=O points toward -c, so that A sits above it
# as in the paper's Fig. 2. Offsets are lattice.rs NN[1]: (sublattice, di, dj, dk).
NN = {
    "1": [(0, 0, 0, 0), (0, 1, -1, 0)], "2": [(0, 0, -1, 0), (0, 1, 0, 0)],
    "A": [(1, 0, 0, 1)], "B": [(1, 0, 0, -1)],
}
CS = 1
center = pos(CS, 0, 0, 0)
CENTER = np.array([center[0] * A, center[1] * A, center[2] * C])
cart = lambda f: np.array([f[0] * A, f[1] * A, f[2] * C])
col = {"1": "#7e22ce", "2": "#dc2626", "A": "#ea580c", "B": "#db2777"}

def molecule(ax, p, s, color, label, proj, size=1.0):
    """A urea molecule as a short bar in its plane with the C=O marked; proj maps 3D -> 2D."""
    # molecular plane: alpha along (1,1,0)/sqrt2 (type-1 direction), beta along (1,-1,0)/sqrt2
    u = np.array([1, 1, 0]) / np.sqrt(2) if s == 0 else np.array([1, -1, 0]) / np.sqrt(2)
    zdir = np.array([0, 0, 1.0]) * (1 if s == 0 else -1)   # C=O direction
    n1, n2 = p - 1.1 * size * u, p + 1.1 * size * u          # the two N
    o = p + 1.25 * size * zdir                                # the O
    for a, b, lw in [(n1, p, 2.5), (n2, p, 2.5), (p, o, 3.5)]:
        ax.plot(*zip(proj(a), proj(b)), color=color, lw=lw, solid_capstyle="round", zorder=3)
    ax.scatter(*proj(o), s=28, color=color, zorder=4)
    if label:
        q = proj(p); c0 = proj(CENTER)
        d = np.array(q) - np.array(c0); d = d / (np.linalg.norm(d) + 1e-9) * 16  # push the label outward
        if label == "B":
            d = np.array([16.0, -4.0])
        ax.annotate(label, q, textcoords="offset points", xytext=tuple(d), fontsize=10, fontweight="bold",
                    color=color, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.0, 3.6))

# left: the molecule's own plane, horizontal axis along [110] (type-1 direction), vertical along c
e_u = np.array([1, -1, 0]) / np.sqrt(2)  # the beta molecule's own plane direction
proj_side = lambda v: (float(np.dot(v, e_u)), float(v[2]))
molecule(ax0, cart(center), CS, "#0f172a", "", proj_side, 1.0)
for t in ["A", "B", "1"]:
    for (s, di, dj, dk) in NN[t]:
        q = cart(pos(s, di, dj, dk))
        molecule(ax0, q, s, col[t], t, proj_side, 1.0)
        ax0.plot(*zip(proj_side(cart(center)), proj_side(q)), color="#94a3b8", lw=0.8, ls=":", zorder=1)
ax0.set_aspect("equal"); ax0.set_xlabel("along [1$\\bar{1}$0] (\u00c5)"); ax0.set_ylabel("along c (\u00c5)")
ax0.set_title("in the molecule's own plane", fontsize=9)
ax0.annotate("O", proj_side(cart(center) + np.array([0, 0, -1.25])), textcoords="offset points", xytext=(6, -3), fontsize=7, color="#0f172a")
ax0.annotate("NH$_2$  NH$_2$", proj_side(cart(center)), textcoords="offset points", xytext=(-22, 6), fontsize=7, color="#0f172a")

# right: view along c, all four in-plane partners (both type-2 are the other orientation)
proj_top = lambda v: (float(v[0]), float(v[1]))
molecule(ax1, cart(center), CS, "#0f172a", "", proj_top, 1.0)
for t in ["1", "2"]:
    for (s, di, dj, dk) in NN[t]:
        q = cart(pos(s, di, dj, dk))
        molecule(ax1, q, s, col[t], t, proj_top, 1.0)
        ax1.plot(*zip(proj_top(cart(center)), proj_top(q)), color="#94a3b8", lw=0.8, ls=":", zorder=1)
ax1.set_aspect("equal"); ax1.set_xlabel("a (\u00c5)"); ax1.set_ylabel("b (\u00c5)")
ax1.set_title("viewed along c", fontsize=9)
for ax in (ax0, ax1):
    ax.margins(0.25)
fig.tight_layout()
fig.savefig(out / "lattice-check.pdf")
print(f"wrote {out}/lattice-check.pdf")
