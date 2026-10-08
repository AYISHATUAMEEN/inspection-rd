#!/usr/bin/env python3
"""02_rd_analysis.py - confirmatory RD analysis exactly as pre-specified in docs/PAP.md (sections 4-6).

Run only after registration (see guard.py), or with --dry-run on a simulated outcome.
Writes output/rd_results[_dryrun].json and .txt.
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats
import rdlib, samples as S, guard

MODE = guard.mode(); TAG = "_dryrun" if MODE == "dry-run" else ""
OUT = os.path.join(S.ROOT, "output"); os.makedirs(OUT, exist_ok=True)
H_PRIMARY, H_ALT = 5.0, (3.0, 8.0)
i = S.load(); h = S.history_pairs(i); p = S.primary(h)
if MODE == "dry-run":
    p["next_score"] = np.where(p.has_next == 1, guard.simulated_outcome(p.inspection_score.values), np.nan)

# fixed-horizon outcome: score at the first listed inspection at least 21 months after the index inspection
hh = h[["ahcp_property_id", "inspection_date", "inspection_score"]].rename(columns={"inspection_date": "d2", "inspection_score": "s2"})
m = p[["ahcp_inspection_id", "ahcp_property_id", "inspection_date"]].merge(hh, on="ahcp_property_id")
m = m[m.d2 >= m.inspection_date + pd.DateOffset(months=21)].sort_values("d2").groupby("ahcp_inspection_id").first()
p["score_21m"] = p.ahcp_inspection_id.map(m.s2)
if MODE == "dry-run": p["score_21m"] = np.where(p.score_21m.notna(), guard.simulated_outcome(p.inspection_score.values, seed=5), np.nan)
p["next_fail"] = np.where(p.next_score.notna(), (np.floor(p.next_score + 0.5) < 60).astype(float), np.nan)
p["index_year"] = p.inspection_date.dt.year.astype(float); p["index_month"] = p.inspection_date.dt.month.astype(float)
x = p.inspection_score.values; cl = p.ahcp_property_id.values
R = {"mode": MODE, "n_index": int(len(p))}; L = [f"RD RESULTS ({MODE})", "=" * 40]

def est(y, c, hbw, d=None, donut=0.0):
    o = rdlib.local_linear_rd(x, y, c, hbw, cluster=cl, d=d, donut=donut)
    o["p"] = float(2 * stats.norm.sf(abs(o["jump"] / o["se"]))); o["ci"] = [o["jump"] - 1.96 * o["se"], o["jump"] + 1.96 * o["se"]]
    return o

# ---- primary family: reduced-form jump in the next score at 80 and 90 (Holm across the two) -----------------
prim = {k: est(p.next_score.values, c, H_PRIMARY, d=p.years_to_next.values) for k, c in S.CUTOFFS.items()}
ps = sorted(prim.items(), key=lambda kv: kv[1]["p"])
holm = {ps[0][0]: min(1.0, 2 * ps[0][1]["p"]), ps[1][0]: min(1.0, max(2 * ps[0][1]["p"], ps[1][1]["p"]))}
for k in prim: prim[k]["p_holm"] = holm[k]
R["primary"] = prim
L.append("\nPrimary: jump in next score at the cutoff (local linear, triangular, h = 5, clustered by property)")
for k, o in prim.items():
    L.append(f"  cutoff {k}: {o['jump']:+.2f} (se {o['se']:.2f}; 95% CI {o['ci'][0]:+.2f} to {o['ci'][1]:+.2f}); p = {o['p']:.4f}; Holm p = {o['p_holm']:.4f}; n = {o['n']:,}")
    L.append(f"     first stage {o['first_stage']:+.3f} years (se {o['first_stage_se']:.3f}); fuzzy estimate per additional year {o['wald']:+.2f} (se {o['wald_se']:.2f})")

# ---- secondary outcomes -------------------------------------------------------------------------------------
R["secondary"] = {}
L.append("\nSecondary outcomes (h = 5)")
for name, y in (("next_fail", p.next_fail.values), ("score_21m", p.score_21m.values), ("has_next", p.has_next.values.astype(float))):
    for k, c in S.CUTOFFS.items():
        o = est(y, c, H_PRIMARY); R["secondary"][f"{name}_{k}"] = o
        L.append(f"  {name}, cutoff {k}: {o['jump']:+.3f} (se {o['se']:.3f}); p = {o['p']:.4f}; n = {o['n']:,}")

# ---- robustness -----------------------------------------------------------------------------------------------
R["robust"] = {}
L.append("\nRobustness of the primary estimates")
for k, c in S.CUTOFFS.items():
    for hbw in H_ALT:
        o = est(p.next_score.values, c, hbw); R["robust"][f"h{hbw:g}_{k}"] = o; L.append(f"  cutoff {k}, h = {hbw:g}: {o['jump']:+.2f} (se {o['se']:.2f})")
    o = rdlib.local_linear_rd(x, p.next_score.values, c, 8.0, cluster=cl, order=2); R["robust"][f"quadratic_h8_{k}"] = o; L.append(f"  cutoff {k}, local quadratic, h = 8: {o['jump']:+.2f} (se {o['se']:.2f})")
    for dn in (0.5, 1.0):
        o = est(p.next_score.values, c, H_PRIMARY, donut=dn); R["robust"][f"donut{dn:g}_{k}"] = o; L.append(f"  cutoff {k}, donut {dn:g}: {o['jump']:+.2f} (se {o['se']:.2f})")
    first = (p.seq_hist.values == 1) | np.isnan(p.prev_score_hist.values)
    o = rdlib.local_linear_rd(x[~first], p.next_score.values[~first], c, H_PRIMARY, cluster=cl[~first]); R["robust"][f"not_first_{k}"] = o
    L.append(f"  cutoff {k}, excluding a property's first listed inspection: {o['jump']:+.2f} (se {o['se']:.2f})")

# ---- attrition and Lee (2009) trimming bounds -------------------------------------------------------------------
# A longer interval leaves more time for a property to leave HUD's files before its next inspection, so the share
# with an observed outcome is lower above each cutoff. Bounds: trim the below-cutoff outcomes (the side with more
# follow-up) by the excess share q, from the top and from the bottom, and re-estimate.
R["lee"] = {}
L.append("\nAttrition and trimming bounds (h = 5)")
for k, c in S.CUTOFFS.items():
    a_ = est(p.has_next.values.astype(float), c, H_PRIMARY); gap = -a_["jump"]; q = max(gap, 0.0) / a_["left_limit"] if a_["left_limit"] > 0 else 0.0
    y = p.next_score.values.copy(); win = (np.abs(x - c) <= H_PRIMARY) & (x < c) & np.isfinite(y)
    lo_cut, hi_cut = np.quantile(y[win], q), np.quantile(y[win], 1 - q)
    y_hi = y.copy(); y_hi[win & (y > hi_cut)] = np.nan      # drop the best below-cutoff outcomes -> largest jump
    y_lo = y.copy(); y_lo[win & (y < lo_cut)] = np.nan      # drop the worst below-cutoff outcomes -> smallest jump
    ub, lb = est(y_hi, c, H_PRIMARY), est(y_lo, c, H_PRIMARY)
    R["lee"][k] = {"attrition_jump": a_["jump"], "attrition_se": a_["se"], "trim_share": q, "lower": lb["jump"], "lower_se": lb["se"], "upper": ub["jump"], "upper_se": ub["se"]}
    L.append(f"  cutoff {k}: attrition jump {a_['jump']:+.3f} (se {a_['se']:.3f}); trim share {q:.3f}; bounds [{lb['jump']:+.2f}, {ub['jump']:+.2f}]; "
             f"95% interval for the bounds [{lb['jump'] - 1.96 * lb['se']:+.2f}, {ub['jump'] + 1.96 * ub['se']:+.2f}]")

# ---- validity: density, predetermined covariates, placebo cutoffs ------------------------------------------------
R["validity"] = {}
L.append("\nValidity checks")
for k, c in S.CUTOFFS.items():
    d = rdlib.density_ratio(x, c, 6.0, 0.5); R["validity"][f"density_{k}"] = d
    L.append(f"  density at {k}: log ratio {d['log_ratio']:+.3f} (se {d['se']:.3f}, z = {d['z']:+.2f})")
    for cov in ("prev_score_hist", "index_year", "index_month", "seq_hist"):
        o = est(p[cov].values.astype(float), c, H_PRIMARY); R["validity"][f"{cov}_{k}"] = o
        L.append(f"  balance, {cov} at {k}: {o['jump']:+.3f} (se {o['se']:.3f}); p = {o['p']:.3f}")
for c in (74.5, 84.5, 94.5):
    o = est(p.next_score.values, c, 4.0); R["validity"][f"placebo_{c}"] = o
    L.append(f"  placebo cutoff {c}: {o['jump']:+.2f} (se {o['se']:.2f}); p = {o['p']:.3f}")

# ---- data-driven bandwidth and robust bias-corrected inference (PAP 6.3), rdrobust / rddensity ---------------
import warnings; warnings.filterwarnings("ignore")
import rdpkgs
R["software"] = rdpkgs.VERSIONS; R["rdrobust"] = {}; R["rddensity"] = {}
L.append("\nrdrobust: MSE-optimal bandwidth, triangular kernel, local linear, clustered by property")
for k, c in S.CUTOFFS.items():
    ok = np.isfinite(p.next_score.values)
    r = rdpkgs.rdrobust(y=p.next_score.values[ok], x=x[ok], c=c, kernel="triangular", p=1, cluster=cl[ok])
    hh = float(r.bws.iloc[0, 0]); conv = float(r.coef.iloc[0, 0]); bc = float(r.coef.iloc[1, 0])
    ci = [float(v) for v in r.ci.iloc[2].values]; pv = float(r.pv.iloc[2, 0])
    R["rdrobust"][k] = {"h": hh, "conventional": conv, "bias_corrected": bc, "robust_ci": ci, "robust_p": pv, "n_eff": [int(v) for v in r.N_h]}
    L.append(f"  cutoff {k}: h = {hh:.2f}; conventional {conv:+.2f}; bias-corrected {bc:+.2f}; robust 95% CI {ci[0]:+.2f} to {ci[1]:+.2f}; robust p = {pv:.4f}")
    d = rdpkgs.rddensity(X=x, c=c)
    tp = float(d.test["p_jk"]); tt = float(d.test["t_jk"])
    R["rddensity"][k] = {"t": tt, "p": tp, "h_left": float(d.h["left"]) if hasattr(d, "h") and hasattr(d.h, "__getitem__") else None}
    L.append(f"  rddensity at {k}: t = {tt:+.2f}, p = {tp:.4f}")

# ---- binned means for the RD figures -------------------------------------------------------------------------------
b = np.floor(x + 0.5); ok = p.next_score.notna().values
R["bins"] = {int(k): [int(v["size"]), float(v["mean"])] for k, v in pd.DataFrame({"b": b[ok], "y": p.next_score.values[ok]}).groupby("b").y.agg(["size", "mean"]).loc[65:99].iterrows()}
json.dump(R, open(os.path.join(OUT, f"rd_results{TAG}.json"), "w"), indent=1, default=float)
open(os.path.join(OUT, f"rd_results{TAG}.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
