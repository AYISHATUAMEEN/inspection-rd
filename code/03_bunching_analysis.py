#!/usr/bin/env python3
"""03_bunching_analysis.py - the score-59 analysis exactly as pre-specified in docs/PAP.md (section 7).

Uses only the distribution of scores (no next-inspection outcome). Gated like the RD script so that the
confirmatory numbers are produced after registration. Writes output/bunching_results[_dryrun].json and .txt.
"""
import json, os
import numpy as np
import rdlib, samples as S, guard

MODE = guard.mode(); TAG = "_dryrun" if MODE == "dry-run" else ""
OUT = os.path.join(S.ROOT, "output"); i = S.load()
if MODE == "dry-run":   # simulated integer scores with a pure relocation of half of 60-70 to 59 in the "NSPIRE" rows
    rng = np.random.default_rng(3); s = np.clip(np.round(100 - rng.gamma(2.0, 6.0, len(i))), 0, 100)
    mv = (i.protocol.values == "NSPIRE") & (s >= 60) & (s <= 70) & (rng.uniform(size=len(i)) < 0.5); s[mv] = 59
    i = i.assign(inspection_score=s)
R = {"mode": MODE}; L = [f"BUNCHING RESULTS ({MODE})", "=" * 40]
def run(lab, g, key):
    b = rdlib.bunching(g.inspection_score.values, 40, 85, 55, 72, degree=5, points=(59,), deficit=(60, 70))
    R[key] = {k: v for k, v in b.items()}
    z = b["excess_minus_missing"] / b["excess_minus_missing_se"] if b["excess_minus_missing_se"] > 0 else float("nan")
    L.append(f"{lab}: n in window {b['n_window']:,}; excess at 59 = {b['excess']:.0f} (se {b['excess_se']:.0f}); missing in 60-70 = {b['missing']:.0f} (se {b['missing_se']:.0f}); "
             f"difference = {b['excess_minus_missing']:.0f} (se {b['excess_minus_missing_se']:.0f}, z = {z:+.2f}); excess ratio {b['excess_ratio']:.2f}")
for prog in ("MF", "PH"):
    run(f"NSPIRE {prog}", i[(i.program == prog) & (i.protocol == "NSPIRE")], f"nspire_{prog}")
    run(f"UPCS 2013-2019 {prog} (no unit-threshold rule; comparison)", i[(i.program == prog) & (i.protocol == "UPCS") & i.inspection_year.between(2013, 2019)], f"upcs_{prog}")
# robustness of the NSPIRE decomposition to polynomial order and excluded window
for prog in ("MF", "PH"):
    g = i[(i.program == prog) & (i.protocol == "NSPIRE")].inspection_score.values
    for deg, (e0, e1) in ((4, (55, 72)), (6, (55, 72)), (5, (53, 74)), (5, (56, 71))):
        b = rdlib.bunching(g, 40, 85, e0, e1, degree=deg, reps=200); R[f"nspire_{prog}_deg{deg}_excl{e0}_{e1}"] = {k: b[k] for k in ("excess", "excess_se", "missing", "missing_se", "excess_minus_missing", "excess_minus_missing_se")}
        L.append(f"NSPIRE {prog}, order {deg}, excluded {e0}-{e1}: excess {b['excess']:.0f}; missing {b['missing']:.0f}; difference {b['excess_minus_missing']:.0f} (se {b['excess_minus_missing_se']:.0f})")
# history sample: density at the passing score (decimal scores)
hs = S.primary(S.history_pairs(S.load())).inspection_score.values
d = rdlib.density_ratio(hs, 59.5, 6.0, 0.5); R["history_density_60"] = d
L.append(f"History sample 2001-2005, density at 59.5: log ratio {d['log_ratio']:+.3f} (se {d['se']:.3f}, z = {d['z']:+.2f})")
# density discontinuities at the three thresholds in integer-score eras: counts just below vs just above
for prog in ("MF", "PH"):
    g = i[(i.program == prog) & (i.protocol == "UPCS") & i.inspection_year.between(2013, 2019)]; sc = np.round(g.inspection_score.values)
    for c in (60, 80, 90):
        lo, hi = int(((sc >= c - 2) & (sc < c)).sum()), int(((sc >= c) & (sc < c + 2)).sum())
        # counterfactual share above from a local linear fit to counts in [c-10, c+9] excluding [c-2, c+1]
        grid = np.arange(c - 10, c + 10); cnt = np.array([(sc == v).sum() for v in grid], float); inc = (grid < c - 2) | (grid > c + 1)
        co = np.polyfit(grid[inc], cnt[inc], 2); cf = np.polyval(co, grid)
        exp_lo, exp_hi = cf[(grid >= c - 2) & (grid < c)].sum(), cf[(grid >= c) & (grid < c + 2)].sum()
        R[f"threshold_{prog}_{c}"] = {"below": lo, "above": hi, "cf_below": float(exp_lo), "cf_above": float(exp_hi)}
        L.append(f"UPCS 2013-2019 {prog}, threshold {c}: observed {lo:,} in [{c-2},{c-1}] and {hi:,} in [{c},{c+1}]; counterfactual {exp_lo:.0f} and {exp_hi:.0f}")
json.dump(R, open(os.path.join(OUT, f"bunching_results{TAG}.json"), "w"), indent=1, default=float)
open(os.path.join(OUT, f"bunching_results{TAG}.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
