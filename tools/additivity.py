"""Bond-additivity test of the Piana-Reyhani-Gale rate tables.

Run: uv run python tools/additivity.py
Fits ln(kd/kc), ln kd and ln kc over the measured cells to per-neighbor terms (1, 2, A, B) plus a constant.
"""
import json, pathlib
import numpy as np

root = pathlib.Path(__file__).resolve().parent.parent
d = json.loads((root / "data/piana_gale_2005_rates.json").read_text())
cols = {"-": (0, 0), "1": (1, 0), "11": (2, 0), "2": (0, 1), "12": (1, 1), "112": (2, 1), "22": (0, 2), "122": (1, 2), "1122": (2, 2)}
rows = {"-": (0, 0), "A": (1, 0), "B": (0, 1), "AB": (1, 1)}
for solv in ["water", "methanol"]:
    X, ys = [], {"ln(kd/kc)": [], "ln kd": [], "ln kc": []}
    for c, (n1, n2) in cols.items():
        for r, (a, b) in rows.items():
            if (c, r) == ("-", "-"):
                continue
            kd, kc = d[f"{solv}_dissolution"][f"{c}/{r}"], d[f"{solv}_crystallization"][f"{c}/{r}"]
            for x, y in zip(kd, kc):
                if x <= 0 or y <= 0 or (x == 1e10 and y == 1e10) or x == y:
                    continue  # unobserved, placeholder, or the "-" normalization
                X.append([n1, n2, a, b, 1.0])
                ys["ln(kd/kc)"].append(np.log(x / y)); ys["ln kd"].append(np.log(x)); ys["ln kc"].append(np.log(y))
    X = np.array(X)
    for lab, y in ys.items():
        y = np.array(y)
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        pred = X @ beta
        r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        rms = np.sqrt(((y - pred) ** 2).mean())
        print(f"{solv:9s} {lab:10s} n={len(y):2d}  R2={r2:.2f}  scatter x{np.exp(rms):.1f}  "
              f"per neighbor (kT): 1={beta[0]:+.2f} 2={beta[1]:+.2f} A={beta[2]:+.2f} B={beta[3]:+.2f}")
