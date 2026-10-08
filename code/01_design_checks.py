#!/usr/bin/env python3
"""01_design_checks.py - feasibility and design diagnostics run BEFORE registration.

By construction this script never compares the outcome (next score) across a cutoff. It reports sample sizes,
the distribution of the running variable, the first stage on the inspection interval, and power computed from
the outcome's residual variance around a smooth global fit. Writes output/design_checks.{txt,json}.
"""
import json, os
import numpy as np, pandas as pd
import rdlib, samples as S

OUT = os.path.join(S.ROOT, "output"); os.makedirs(OUT, exist_ok=True)
L, J = [], {}
def w(*a): L.append(" ".join(str(x) for x in a))

i = S.load(); h = S.history_pairs(i); p = S.primary(h)
w("DESIGN CHECKS (pre-registration; no outcome contrasts)\n" + "=" * 54)
w(f"History-file multifamily inspections: {len(h):,}; properties: {h.ahcp_property_id.nunique():,}")
w(f"Index inspections {S.RULE_START.date()} to {S.INDEX_END.date()}: {len(p):,}; properties: {p.ahcp_property_id.nunique():,}")
w(f"  with a next inspection listed: {p.has_next.sum():,} ({100 * p.has_next.mean():.1f}%)")
J["n_index"] = int(len(p)); J["n_index_props"] = int(p.ahcp_property_id.nunique()); J["share_has_next"] = float(p.has_next.mean())
x = p.inspection_score.values
w(f"  score decimals: {100 * np.mean(np.abs(x - np.round(x)) > 1e-9):.1f}% of index scores are non-integer (continuous running variable)")
J["share_noninteger"] = float(np.mean(np.abs(x - np.round(x)) > 1e-9))

w("\n1. Observations near each cutoff (index inspections)")
for name, c in list(S.CUTOFFS.items()) + [("60", 59.5)]:
    row = {}
    for hbw in (2, 3, 5, 8, 10):
        m = np.abs(x - c) <= hbw; row[hbw] = [int(((x < c) & m).sum()), int(((x >= c) & m).sum())]
    w(f"   cutoff {name} (c = {c}): " + "; ".join(f"h={k}: {v[0]:,} below / {v[1]:,} above" for k, v in row.items()))
    J[f"n_near_{name}"] = row

w("\n2. Running-variable density at each cutoff (log ratio above/below; local linear on 0.5-point bins, h = 6)")
for name, c in list(S.CUTOFFS.items()) + [("60", 59.5)]:
    d = rdlib.density_ratio(x, c, 6.0, 0.5); J[f"density_{name}"] = d
    w(f"   cutoff {name}: log ratio {d['log_ratio']:+.3f} (se {d['se']:.3f}, z = {d['z']:+.2f}); implied excess above = {100 * (np.exp(d['log_ratio']) - 1):+.1f}%")
w("   Counts in 0.5-point bins, 77.0 to 82.0 and 87.0 to 92.0:")
for lo in (77.0, 87.0):
    e = np.arange(lo, lo + 5.01, 0.5); cnt = np.histogram(x, bins=e)[0]
    w("     " + "  ".join(f"[{a:.1f}) {n}" for a, n in zip(e[:-1], cnt)))

w("\n3. First stage: years to the next listed inspection (treatment), by rounded index score")
pn = p[p.has_next == 1]; r = np.floor(pn.inspection_score + 0.5).astype(int)
tab = pn.groupby(r).years_to_next.agg(["size", "mean", "median"])
for s_ in range(70, 100): w(f"   score {s_}: n={int(tab.loc[s_, 'size']):>5}  mean {tab.loc[s_, 'mean']:.2f}  median {tab.loc[s_, 'median']:.2f}")
J["first_stage_by_score"] = {int(k): [int(v["size"]), float(v["mean"]), float(v["median"])] for k, v in tab.loc[60:99].iterrows()}
w("   Local-linear jump in years to next inspection (triangular kernel, clustered by property):")
for name, c in S.CUTOFFS.items():
    for hbw in (3, 5, 8):
        o = rdlib.local_linear_rd(pn.inspection_score.values, pn.years_to_next.values, c, hbw, cluster=pn.ahcp_property_id.values)
        w(f"     cutoff {name}, h={hbw}: {o['jump']:+.3f} years (se {o['se']:.3f}); n = {o['n']:,}; mean just below {o['left_limit']:.2f}")
        J[f"first_stage_{name}_h{hbw}"] = o
    for lab, lo, hi in (("on-schedule share below (interval <= 1.5y)" if name == "80" else "on-schedule share below (interval <= 2.5y)", None, None),):
        thr = 1.5 if name == "80" else 2.5
        o = rdlib.local_linear_rd(pn.inspection_score.values, (pn.years_to_next.values > thr).astype(float), c, 5, cluster=pn.ahcp_property_id.values)
        w(f"     cutoff {name}, h=5: jump in P(interval > {thr} y) = {o['jump']:+.3f} (se {o['se']:.3f}); just below {o['left_limit']:.3f}")
        J[f"first_stage_prob_{name}"] = o

w("\n4. Power (outcome variance only; no outcome contrast across the cutoff)")
rng = np.random.default_rng(7)
ok = pn.next_score.notna()
xs, ys = pn.inspection_score.values[ok], pn.next_score.values[ok]
m = (xs >= 60) & (xs <= 100)
coef = np.polyfit(xs[m], ys[m], 4); resid = ys[m] - np.polyval(coef, xs[m])   # smooth global fit, no cutoff terms
sd = float(resid.std()); w(f"   residual SD of next score around a quartic in the index score: {sd:.2f} points")
J["resid_sd"] = sd
for name, c in S.CUTOFFS.items():
    for hbw in (3, 5, 8):
        ses = []
        for _ in range(200):
            ysim = np.polyval(coef, pn.inspection_score.values) + rng.normal(0, sd, len(pn))
            ses.append(rdlib.local_linear_rd(pn.inspection_score.values, ysim, c, hbw, cluster=pn.ahcp_property_id.values)["jump"])
        se = float(np.std(ses)); w(f"   cutoff {name}, h={hbw}: simulated SE of the reduced-form jump {se:.2f}; minimum detectable effect (80% power, 5% two-sided) {2.8 * se:.2f} points")
        J[f"mde_{name}_h{hbw}"] = 2.8 * se

w("\n4b. Attrition: jump in P(next inspection listed) at each cutoff, by last index date (h = 5)")
for end in ("2005-12-31", "2005-06-30", "2004-12-31", "2004-06-30", "2003-12-31"):
    q = p[p.inspection_date <= end]
    row = []
    for name, c in S.CUTOFFS.items():
        o = rdlib.local_linear_rd(q.inspection_score.values, q.has_next.values.astype(float), c, 5.0, cluster=q.ahcp_property_id.values)
        row.append(f"cutoff {name}: {o['jump']:+.3f} (se {o['se']:.3f}), just below {o['left_limit']:.3f}"); J[f"attrition_{name}_{end}"] = o
    w(f"   index <= {end} (n = {len(q):,}): " + "; ".join(row))
w("   Among index inspections without a next listing, share whose property is absent from all later HUD files:")
later = set(i[(i.program == "MF") & (i.last_vintage > 2011)].ahcp_property_id)
nn = p[p.has_next == 0]; w(f"     {100 * (~nn.ahcp_property_id.isin(later)).mean():.1f}% of {len(nn):,}")
J["attrition_absent_later_share"] = float((~nn.ahcp_property_id.isin(later)).mean())

w("\n4c. Predetermined covariates: jump at each cutoff (h = 5)")
p["index_year"] = p.inspection_date.dt.year.astype(float); p["index_month"] = p.inspection_date.dt.month.astype(float)
for cov in ("prev_score_hist", "seq_hist", "index_year", "index_month"):
    row = []
    for name, c in S.CUTOFFS.items():
        o = rdlib.local_linear_rd(p.inspection_score.values, p[cov].values.astype(float), c, 5.0, cluster=p.ahcp_property_id.values)
        row.append(f"cutoff {name}: {o['jump']:+.3f} (se {o['se']:.3f})"); J[f"balance_{cov}_{name}"] = o
    w(f"   {cov}: " + "; ".join(row))

w("\n5. Later-era samples (integer scores; next inspection only partly observed)")
for prog in ("MF", "PH"):
    for prot, lab in (("UPCS", "UPCS 2013-2019"), ("NSPIRE", "NSPIRE")):
        g = i[(i.program == prog) & (i.protocol == prot) & (i.inspection_year >= 2013)]
        if prot == "UPCS": g = g[g.inspection_year <= 2019]
        sc = g.inspection_score.round()
        w(f"   {prog} {lab}: n={len(g):,}; at 59: {(sc == 59).sum():,}; at 60: {(sc == 60).sum():,}; 79: {(sc == 79).sum():,}; 80: {(sc == 80).sum():,}; 89: {(sc == 89).sum():,}; 90: {(sc == 90).sum():,}")
        J[f"late_{prog}_{prot}"] = {"n": int(len(g)), **{str(k): int((sc == k).sum()) for k in (58, 59, 60, 61, 79, 80, 89, 90)}}
w(f"\nLatest inspection date in AHCP v1.0: {i.inspection_date.max().date()} (no inspections on or after 2026-10-01 are available yet)")
J["latest_date"] = str(i.inspection_date.max().date())
open(os.path.join(OUT, "design_checks.txt"), "w").write("\n".join(L) + "\n"); json.dump(J, open(os.path.join(OUT, "design_checks.json"), "w"), indent=1)
print("\n".join(L))
