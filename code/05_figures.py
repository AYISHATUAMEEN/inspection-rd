#!/usr/bin/env python3
"""05_figures.py - figures for the paper (paper/figs/*.pdf) from the data and output/*.json."""
import json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rdlib, samples as S

F = os.path.join(S.ROOT, "paper", "figs"); os.makedirs(F, exist_ok=True)
RD = json.load(open(os.path.join(S.ROOT, "output", "rd_results.json"))); BU = json.load(open(os.path.join(S.ROOT, "output", "bunching_results.json")))
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#9a9993",
                     "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "legend.fontsize": 8, "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.03})
p = S.primary(S.history_pairs(S.load())); pn = p[p.has_next == 1]
def grid(ax): ax.yaxis.grid(True, color=GRID, linewidth=0.5); ax.set_axisbelow(True)

def rdpanel(ax, y, label, ylab, cuts=(79.5, 89.5), lo=70, hi=99.5, h=5.0, binw=0.5):
    x = pn.inspection_score.values; m = (x >= lo) & (x < hi); e = np.arange(lo, hi + 0.01, binw)
    idx = np.digitize(x[m], e) - 1; df = pd.DataFrame({"b": idx, "y": y[m]}).groupby("b").y.agg(["mean", "size"])
    ax.scatter(e[df.index] + binw / 2, df["mean"], s=np.clip(df["size"] / 25, 3, 18), color=BLUE, alpha=0.75, linewidths=0, label="Bin mean (0.5 points)")
    for c in cuts:
        for side in (-1, 1):
            zz = np.linspace(0, h, 30) * side; xs = c + zz
            sel = (np.abs(x - c) <= h) & ((x >= c) if side > 0 else (x < c)); w = 1 - np.abs(x[sel] - c) / h
            b = np.polyfit(x[sel] - c, y[sel], 1, w=np.sqrt(w)); ax.plot(xs, np.polyval(b, zz), color=ORANGE, linewidth=1.6)
        ax.axvline(c, color=MUTED, linestyle=(0, (3, 3)), linewidth=0.7)
    ax.set_xlabel("Index inspection score"); ax.set_ylabel(ylab); ax.set_title(label); grid(ax)

# Fig 1 first stage
fig, ax = plt.subplots(figsize=(6.3, 2.9)); rdpanel(ax, pn.years_to_next.values, "", "Years to next inspection")
ax.text(79.7, 3.6, "two-year\nschedule", fontsize=7.5, color=MUTED); ax.text(89.7, 4.0, "three-year\nschedule", fontsize=7.5, color=MUTED); ax.text(72, 2.2, "annual\nschedule", fontsize=7.5, color=MUTED)
fig.savefig(os.path.join(F, "first_stage.pdf")); plt.close(fig)

# Fig 2 reduced form at both cutoffs
fig, ax = plt.subplots(figsize=(6.3, 2.9)); rdpanel(ax, pn.next_score.values, "", "Score at next inspection")
ax.legend(frameon=False, loc="lower right"); fig.savefig(os.path.join(F, "reduced_form.pdf")); plt.close(fig)

# Fig 3 density of index score
fig, ax = plt.subplots(figsize=(6.3, 2.4)); x = p.inspection_score.values; e = np.arange(50, 100.01, 0.5)
ax.hist(x, bins=e, color=BLUE, linewidth=0); [ax.axvline(c, color=MUTED, linestyle=(0, (3, 3)), linewidth=0.7) for c in (59.5, 79.5, 89.5)]
ax.set_xlabel("Index inspection score (0.5-point bins)"); ax.set_ylabel("Index inspections"); grid(ax); fig.savefig(os.path.join(F, "density.pdf")); plt.close(fig)

# Fig 4 estimates across specifications
specs = [("h = 5 (primary)", "primary", None), ("h = 3", "robust", "h3"), ("h = 8", "robust", "h8"), ("Quadratic, h = 8", "robust", "quadratic_h8"),
         ("Donut 0.5", "robust", "donut0.5"), ("Donut 1", "robust", "donut1"), ("Excl. first listed", "robust", "not_first")]
fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.6), sharey=True)
for ax, (k, title) in zip(axes, (("80", "(a) 80 cutoff"), ("90", "(b) 90 cutoff"))):
    rows = []
    for lab, grp, key in specs:
        o = RD["primary"][k] if grp == "primary" else RD["robust"][f"{key}_{k}"]; rows.append((lab, o["jump"], o["se"]))
    rr = RD["rdrobust"][k]; rows.append(("rdrobust (MSE-opt. h)", rr["bias_corrected"], (rr["robust_ci"][1] - rr["robust_ci"][0]) / 3.92))
    lb = RD["lee"][k]; yy = np.arange(len(rows))[::-1]
    for (lab, est, se), yv in zip(rows, yy):
        ax.errorbar(est, yv, xerr=1.96 * se, fmt="o", color=ORANGE if lab.startswith("h = 5") else BLUE, ms=3.5, elinewidth=1, capsize=0)
    ax.plot([lb["lower"], lb["upper"]], [-1, -1], color=INK, linewidth=2.2); rows.append(("Trimming bounds", None, None))
    ax.axvline(0, color=MUTED, linestyle=(0, (3, 3)), linewidth=0.7); ax.set_title(title); ax.set_xlabel("Jump in next score (points)"); grid(ax)
    ax.set_yticks(list(yy) + [-1]); ax.set_yticklabels([r[0] for r in rows[:-1]] + ["Trimming bounds"])
fig.tight_layout(w_pad=1.2); fig.savefig(os.path.join(F, "robustness.pdf")); plt.close(fig)

# Fig 5 bunching: NSPIRE observed vs counterfactual, MF and PH
fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.5))
for ax, key, title in ((axes[0], "nspire_MF", "(a) Multifamily, NSPIRE"), (axes[1], "nspire_PH", "(b) Public housing, NSPIRE")):
    b = BU[key]; g = np.array(b["grid"]); ax.bar(g, b["observed"], width=0.8, color=BLUE, linewidth=0, label="Observed")
    ax.plot(g, b["counterfactual"], color=ORANGE, linewidth=1.5, label="Counterfactual"); ax.axvspan(54.5, 72.5, color=GRID, alpha=0.5, linewidth=0)
    ax.set_title(title); ax.set_xlabel("Inspection score"); grid(ax)
axes[0].set_ylabel("Inspections"); axes[0].legend(frameon=False, loc="upper left"); fig.tight_layout(w_pad=1.0); fig.savefig(os.path.join(F, "bunching.pdf")); plt.close(fig)
print("figures written")
