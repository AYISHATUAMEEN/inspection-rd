#!/usr/bin/env python3
"""06_paper_numbers.py - LaTeX macros and tables for the paper, from output/*.json. No number is typed by hand."""
import json, os, re
import numpy as np, pandas as pd
from scipy import stats
import samples as S
O = os.path.join(S.ROOT, "output"); P = os.path.join(S.ROOT, "paper"); T = os.path.join(P, "tables"); os.makedirs(T, exist_ok=True)
RD = json.load(open(os.path.join(O, "rd_results.json"))); BU = json.load(open(os.path.join(O, "bunching_results.json"))); DC = json.load(open(os.path.join(O, "design_checks.json")))
N = {}
def neg(s): return "$-$" + s[1:] if s.startswith("-") else s
f = lambda v, d=2: neg(f"{v:.{d}f}"); fs = lambda v, d=2: neg(f"{v:+.{d}f}") if v < 0 else f"+{v:.{d}f}"
c = lambda v: neg(f"{int(round(v)):,}")
def put(k, v): assert re.fullmatch(r"[A-Za-z]+", k), k; N[k] = v
W = {"80": "Eighty", "90": "Ninety"}
for k in ("80", "90"):
    o = RD["primary"][k]; w = W[k]
    put(f"jump{w}", fs(o["jump"])); put(f"se{w}", f(o["se"])); put(f"ciLo{w}", f(o["ci"][0])); put(f"ciHi{w}", f(o["ci"][1])); put(f"p{w}", f(o["p"], 3)); put(f"pHolm{w}", f(o["p_holm"], 3))
    put(f"n{w}", c(o["n"])); put(f"fs{w}", f(o["first_stage"])); put(f"fsSe{w}", f(o["first_stage_se"], 3)); put(f"wald{w}", fs(o["wald"])); put(f"waldSe{w}", f(o["wald_se"]))
    put(f"leftLimit{w}", f(o["left_limit"], 1)); put(f"ciLoSd{w}", f(abs(o["ci"][0]) / DC["resid_sd"]))
    for sec in ("next_fail", "score_21m", "has_next"):
        s_ = RD["secondary"][f"{sec}_{k}"]; key = {"next_fail": "Fail", "score_21m": "Fixed", "has_next": "Attr"}[sec]
        mult = 100 if sec in ("next_fail", "has_next") else 1
        put(f"sec{key}{w}", fs(s_["jump"] * mult, 1 if mult == 100 else 2)); put(f"sec{key}Se{w}", f(s_["se"] * mult, 1 if mult == 100 else 2)); put(f"sec{key}P{w}", f(s_["p"], 3))
    for r in ("h3", "h8", "quadratic_h8", "donut0.5", "donut1", "not_first"):
        o2 = RD["robust"][f"{r}_{k}"]; nm = {"h3": "HThree", "h8": "HEight", "quadratic_h8": "Quad", "donut0.5": "DonutHalf", "donut1": "DonutOne", "not_first": "NotFirst"}[r]
        put(f"rob{nm}{w}", fs(o2["jump"])); put(f"rob{nm}Se{w}", f(o2["se"])); put(f"rob{nm}P{w}", f(2 * stats.norm.sf(abs(o2["jump"] / o2["se"])), 3))
        put(f"rob{nm}CiLo{w}", f(o2["jump"] - 1.96 * o2["se"])); put(f"rob{nm}CiHi{w}", f(o2["jump"] + 1.96 * o2["se"]))
    rr = RD["rdrobust"][k]; put(f"rdrH{w}", f(rr["h"])); put(f"rdrConv{w}", fs(rr["conventional"])); put(f"rdrBc{w}", fs(rr["bias_corrected"])); put(f"rdrLo{w}", f(rr["robust_ci"][0])); put(f"rdrHi{w}", f(rr["robust_ci"][1])); put(f"rdrP{w}", f(rr["robust_p"], 3))
    d = RD["rddensity"][k]; put(f"rddT{w}", fs(d["t"])); put(f"rddP{w}", f(d["p"], 3))
    dv = RD["validity"][f"density_{k}"]; put(f"densZ{w}", fs(dv["z"])); put(f"densLr{w}", fs(dv["log_ratio"], 3))
    lb = RD["lee"][k]; put(f"leeLo{w}", fs(lb["lower"])); put(f"leeHi{w}", fs(lb["upper"])); put(f"leeCiLo{w}", f(lb["lower"] - 1.96 * lb["lower_se"])); put(f"leeCiHi{w}", f(lb["upper"] + 1.96 * lb["upper_se"])); put(f"trim{w}", f(100 * lb["trim_share"], 1))
    for cov, nm in (("prev_score_hist", "Prev"), ("index_year", "Year"), ("index_month", "Month"), ("seq_hist", "Seq")):
        o3 = RD["validity"][f"{cov}_{k}"]; put(f"bal{nm}{w}", fs(o3["jump"])); put(f"bal{nm}Se{w}", f(o3["se"])); put(f"bal{nm}P{w}", f(o3["p"], 2))
for cpl, nm in (("74.5", "SeventyFour"), ("84.5", "EightyFour"), ("94.5", "NinetyFour")):
    o = RD["validity"][f"placebo_{cpl}"]; put(f"plac{nm}", fs(o["jump"])); put(f"placSe{nm}", f(o["se"])); put(f"placP{nm}", f(o["p"], 2))
put("nIndex", c(DC["n_index"])); put("nIndexProps", c(DC["n_index_props"])); put("shareNext", f(100 * DC["share_has_next"], 1)); put("residSd", f(DC["resid_sd"], 1))
put("mdeEighty", f(DC["mde_80_h5"], 1)); put("mdeNinety", f(DC["mde_90_h5"], 1)); put("absentLater", f(100 * DC["attrition_absent_later_share"], 0))
put("histDensZ", fs(BU["history_density_60"]["z"])); put("histDensLr", fs(BU["history_density_60"]["log_ratio"], 3)); put("histDensPct", f(100 * (np.exp(BU["history_density_60"]["log_ratio"]) - 1), 0))
# H2
for prog, nm in (("MF", "Mf"), ("PH", "Ph")):
    b = BU[f"nspire_{prog}"]; z = b["excess_minus_missing"] / b["excess_minus_missing_se"]; pv = 2 * stats.norm.sf(abs(z))
    put(f"exc{nm}", c(b["excess"])); put(f"excSe{nm}", c(b["excess_se"])); put(f"mis{nm}", c(b["missing"])); put(f"misSe{nm}", c(b["missing_se"]))
    put(f"dif{nm}", c(b["excess_minus_missing"])); put(f"difSe{nm}", c(b["excess_minus_missing_se"])); put(f"difZ{nm}", fs(z)); put(f"difP{nm}", f(pv, 2)); put(f"excRatio{nm}", f(b["excess_ratio"], 1)); N[f"_p_{prog}"] = pv
    u = BU[f"upcs_{prog}"]; put(f"upcsExc{nm}", c(u["excess"])); put(f"upcsExcSe{nm}", c(u["excess_se"])); put(f"upcsDif{nm}", c(u["excess_minus_missing"])); put(f"upcsDifSe{nm}", c(u["excess_minus_missing_se"])); put(f"upcsMis{nm}", c(u["missing"]))
ps = sorted([("MF", N.pop("_p_MF")), ("PH", N.pop("_p_PH"))], key=lambda t: t[1]); hp = {ps[0][0]: min(1, 2 * ps[0][1])}; hp[ps[1][0]] = min(1, max(hp[ps[0][0]], ps[1][1]))
put("difPHolmMf", f(hp["MF"], 2)); put("difPHolmPh", f(hp["PH"], 2))
rows = []
for prog, lab in (("MF", "Multifamily"), ("PH", "Public housing")):
    for deg, e0, e1 in ((5, 55, 72), (4, 55, 72), (6, 55, 72), (5, 53, 74), (5, 56, 71)):
        b = BU[f"nspire_{prog}"] if (deg, e0, e1) == (5, 55, 72) else BU[f"nspire_{prog}_deg{deg}_excl{e0}_{e1}"]
        rows.append([lab if deg == 5 and e0 == 55 else "", deg, f"{e0}--{e1}", f"{c(b['excess'])} ({c(b['excess_se'])})", f"{c(b['missing'])} ({c(b['missing_se'])})", f"{c(b['excess_minus_missing'])} ({c(b['excess_minus_missing_se'])})"])
    u = BU[f"upcs_{prog}"]; rows.append([f"\\quad UPCS 2013--19 (comparison)", 5, "55--72", f"{c(u['excess'])} ({c(u['excess_se'])})", f"{c(u['missing'])} ({c(u['missing_se'])})", f"{c(u['excess_minus_missing'])} ({c(u['excess_minus_missing_se'])})"])
    if prog == "MF": rows.append("\\midrule\n")
with open(os.path.join(T, "bunching.tex"), "w") as fh:
    fh.write("\\begin{tabular}{llcrrr}\n\\toprule\nSample & Order & Excluded & Excess at 59 & Missing in 60--70 & Difference \\\\\n\\midrule\n")
    for r in rows: fh.write(r if isinstance(r, str) else " & ".join(str(x) for x in r) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
rows = []
for prog, lab in (("MF", "Multifamily"), ("PH", "Public housing")):
    for t_ in (60, 80, 90):
        r = BU[f"threshold_{prog}_{t_}"]; rows.append([lab if t_ == 60 else "", t_, c(r["below"]), c(r["cf_below"]), f(100 * (r["below"] / r["cf_below"] - 1), 0), c(r["above"]), c(r["cf_above"]), f(100 * (r["above"] / r["cf_above"] - 1), 0)])
        put(f"thr{'Mf' if prog == 'MF' else 'Ph'}{ {60: 'Sixty', 80: 'Eighty', 90: 'Ninety'}[t_]}Below", f(100 * (r["below"] / r["cf_below"] - 1), 0)); put(f"thr{'Mf' if prog == 'MF' else 'Ph'}{ {60: 'Sixty', 80: 'Eighty', 90: 'Ninety'}[t_]}Above", f(100 * (r["above"] / r["cf_above"] - 1), 0))
with open(os.path.join(T, "thresholds.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrrrr}\n\\toprule\n& & \\multicolumn{3}{c}{Two scores below} & \\multicolumn{3}{c}{Two scores at or above} \\\\\n\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\nProgram & Threshold & Observed & Counterf. & Diff. (\\%) & Observed & Counterf. & Diff. (\\%) \\\\\n\\midrule\n")
    for r in rows: fh.write(" & ".join(str(x) for x in r) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
# main RD table
def cell(o): return f"{fs(o['jump'])}", f"({f(o['se'])})"
lines = []
for lab, get in (("Primary: local linear, $h=5$", lambda k: RD["primary"][k]), ("Bandwidth $h=3$", lambda k: RD["robust"][f"h3_{k}"]), ("Bandwidth $h=8$", lambda k: RD["robust"][f"h8_{k}"]),
                 ("Local quadratic, $h=8$", lambda k: RD["robust"][f"quadratic_h8_{k}"]), ("Donut, $|x-c|\\geq0.5$", lambda k: RD["robust"][f"donut0.5_{k}"]), ("Donut, $|x-c|\\geq1$", lambda k: RD["robust"][f"donut1_{k}"]),
                 ("Excluding first listed inspection", lambda k: RD["robust"][f"not_first_{k}"])):
    a, b = cell(get("80")), cell(get("90")); lines.append(f"{lab} & {a[0]} & {b[0]} \\\\\n & {a[1]} & {b[1]} \\\\\n")
r8, r9 = RD["rdrobust"]["80"], RD["rdrobust"]["90"]
lines.append(f"\\texttt{{rdrobust}}, bias-corrected & {fs(r8['bias_corrected'])} & {fs(r9['bias_corrected'])} \\\\\n\\quad robust 95\\% CI; MSE-optimal $h$ & [{f(r8['robust_ci'][0])}, {f(r8['robust_ci'][1])}]; {f(r8['h'])} & [{f(r9['robust_ci'][0])}, {f(r9['robust_ci'][1])}]; {f(r9['h'])} \\\\\n")
l8, l9 = RD["lee"]["80"], RD["lee"]["90"]
lines.append(f"Trimming bounds (attrition) & [{fs(l8['lower'])}, {fs(l8['upper'])}] & [{fs(l9['lower'])}, {fs(l9['upper'])}] \\\\\n")
with open(os.path.join(T, "rd_main.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lcc}\n\\toprule\n& 80 cutoff & 90 cutoff \\\\\n\\midrule\n" + "".join(lines) + "\\midrule\n")
    fh.write(f"Observations, primary & {c(RD['primary']['80']['n'])} & {c(RD['primary']['90']['n'])} \\\\\nFirst stage (years), primary & {f(RD['primary']['80']['first_stage'])} ({f(RD['primary']['80']['first_stage_se'], 3)}) & {f(RD['primary']['90']['first_stage'])} ({f(RD['primary']['90']['first_stage_se'], 3)}) \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
# summary statistics of primary sample
p = S.primary(S.history_pairs(S.load())); pn = p[p.has_next == 1]
rows = [("Index score", p.inspection_score), ("Next score", pn.next_score), ("Years to next inspection", pn.years_to_next), ("Previous score (where listed)", p.prev_score_hist.dropna()), ("Inspection year", p.inspection_date.dt.year.astype(float))]
with open(os.path.join(T, "summary.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrr}\n\\toprule\nVariable & $n$ & Mean & SD & P10 & P90 \\\\\n\\midrule\n")
    for lab, s in rows: fh.write(f"{lab} & {c(len(s))} & {f(s.mean())} & {f(s.std())} & {f(s.quantile(.1))} & {f(s.quantile(.9))} \\\\\n")
    fh.write(f"Share below 60 at index & {c(len(p))} & {f(100*(p.inspection_score < 59.5).mean(),1)}\\% & & & \\\\\nShare 80 or above at index & {c(len(p))} & {f(100*(p.inspection_score >= 79.5).mean(),1)}\\% & & & \\\\\nShare 90 or above at index & {c(len(p))} & {f(100*(p.inspection_score >= 89.5).mean(),1)}\\% & & & \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
# validity table: balance, placebo, density
rows = []
for cov, lab in (("prev_score_hist", "Previous score"), ("seq_hist", "Inspection sequence number"), ("index_year", "Index year"), ("index_month", "Index month")):
    a8, a9 = RD["validity"][f"{cov}_80"], RD["validity"][f"{cov}_90"]
    rows.append([lab, f"{fs(a8['jump'])} ({f(a8['se'])})", f(a8["p"], 2), f"{fs(a9['jump'])} ({f(a9['se'])})", f(a9["p"], 2)])
a8, a9 = RD["secondary"]["has_next_80"], RD["secondary"]["has_next_90"]
rows.append(["Next inspection listed (attrition)", f"{fs(a8['jump'], 3)} ({f(a8['se'], 3)})", f(a8["p"], 3), f"{fs(a9['jump'], 3)} ({f(a9['se'], 3)})", f(a9["p"], 3)])
d8, d9 = RD["validity"]["density_80"], RD["validity"]["density_90"]
rows.append(["Density, binned log ratio", f"{fs(d8['log_ratio'], 3)} ({f(d8['se'], 3)})", f(2 * stats.norm.sf(abs(d8['z'])), 3), f"{fs(d9['log_ratio'], 3)} ({f(d9['se'], 3)})", f(2 * stats.norm.sf(abs(d9['z'])), 3)])
r8, r9 = RD["rddensity"]["80"], RD["rddensity"]["90"]
rows.append(["Density, \\texttt{rddensity} $t$", fs(r8["t"]), f(r8["p"], 3), fs(r9["t"]), f(r9["p"], 3)])
with open(os.path.join(T, "validity.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrr}\n\\toprule\n& \\multicolumn{2}{c}{80 cutoff} & \\multicolumn{2}{c}{90 cutoff} \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\nOutcome & Jump (SE) & $p$ & Jump (SE) & $p$ \\\\\n\\midrule\n")
    for r in rows: fh.write(" & ".join(r) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
rows = []
for cpl in ("74.5", "84.5", "94.5"):
    o = RD["validity"][f"placebo_{cpl}"]; rows.append([cpl, f"{fs(o['jump'])} ({f(o['se'])})", f(o["p"], 2), c(o["n"])])
with open(os.path.join(T, "placebo.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrr}\n\\toprule\nPlacebo cutoff & Jump (SE) & $p$ & $n$ \\\\\n\\midrule\n")
    for r in rows: fh.write(" & ".join(r) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
rows = []
for sec, lab, mult, dd in (("next_fail", "Next inspection fails (pp)", 100, 1), ("score_21m", "Fixed-horizon score (points)", 1, 2)):
    a8, a9 = RD["secondary"][f"{sec}_80"], RD["secondary"][f"{sec}_90"]
    rows.append([lab, f"{fs(a8['jump'] * mult, dd)} ({f(a8['se'] * mult, dd)})", f(a8["p"], 3), c(a8["n"]), f"{fs(a9['jump'] * mult, dd)} ({f(a9['se'] * mult, dd)})", f(a9["p"], 3), c(a9["n"])])
for k, lab in (("80", "Fuzzy: per additional year (points)"),):
    pass
a8, a9 = RD["primary"]["80"], RD["primary"]["90"]
rows.append(["Fuzzy: points per additional year", f"{fs(a8['wald'])} ({f(a8['wald_se'])})", f(2 * stats.norm.sf(abs(a8['wald'] / a8['wald_se'])), 3), c(a8["n"]), f"{fs(a9['wald'])} ({f(a9['wald_se'])})", f(2 * stats.norm.sf(abs(a9['wald'] / a9['wald_se'])), 3), c(a9["n"])])
with open(os.path.join(T, "secondary.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n& \\multicolumn{3}{c}{80 cutoff} & \\multicolumn{3}{c}{90 cutoff} \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\nOutcome & Jump (SE) & $p$ & $n$ & Jump (SE) & $p$ & $n$ \\\\\n\\midrule\n")
    for r in rows: fh.write(" & ".join(r) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
SIM = json.load(open(os.path.join(O, "simulation_tests.json")))
with open(os.path.join(T, "simulation.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lll}\n\\toprule\nCheck & Result & Criterion met \\\\\n\\midrule\n")
    for r in SIM: fh.write(f"{r['test'].capitalize().replace(' rd ', ' RD ')} & {r['detail'].replace('-', '$-$')} & {'yes' if r['pass'] else 'no'} \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
put("simAllPass", "all" if all(r["pass"] for r in SIM) else "not all")
with open(os.path.join(P, "numbers.tex"), "w") as fh:
    for k, v in sorted(N.items()): fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
json.dump(N, open(os.path.join(P, "paper_numbers.json"), "w"), indent=1); print(len(N), "macros")
