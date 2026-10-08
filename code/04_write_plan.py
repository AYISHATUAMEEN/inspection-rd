#!/usr/bin/env python3
"""04_write_plan.py - write docs/PAP.md from output/design_checks.json so that every number in the plan is computed."""
import json, os, datetime
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = json.load(open(os.path.join(ROOT, "output", "design_checks.json"))); A = json.load(open(os.path.join(ROOT, "AUTHORS.json")))["authors"][0]
TODAY = datetime.date.today().isoformat()
n = lambda v: f"{v:,}"
fs80, fs90 = J["first_stage_80_h5"], J["first_stage_90_h5"]; pr80, pr90 = J["first_stage_prob_80"], J["first_stage_prob_90"]
d80, d90, d60 = J["density_80"], J["density_90"], J["density_60"]
at80, at90 = J["attrition_80_2005-12-31"], J["attrition_90_2005-12-31"]
near = lambda k, h: f"{n(J['near'][k][h][0])} below and {n(J['near'][k][h][1])} above" if False else f"{n(J[f'n_near_{k}'][str(h)][0])} below and {n(J[f'n_near_{k}'][str(h)][1])} above"
bal = lambda cov, k: f"{J[f'balance_{cov}_{k}']['jump']:+.2f} (se {J[f'balance_{cov}_{k}']['se']:.2f})"
lm = J["late_MF_NSPIRE"]; lp = J["late_PH_NSPIRE"]; um = J["late_MF_UPCS"]

doc = f"""# Pre-Analysis Plan

**Do Inspection Thresholds Improve Housing Conditions? Regression Discontinuity Evidence from HUD's Physical Inspection Schedule**

{A['name']} ({A['affiliation']}; ORCID {A['orcid']})

Plan version 1.0, {TODAY}. Template: OSF Secondary Data Preregistration. This plan is written before any outcome contrast has been estimated. Section 3 states exactly what the author has and has not seen.

## 1. Study information

### 1.1 Background

HUD inspects assisted multifamily properties on a schedule set by the previous inspection score. Under 24 CFR 200.857(b), in force from 8 January 2001 (65 FR 77230), a property scoring 90 or above is inspected once every three years, a property scoring 80 to 89 every two years, and a property scoring below 80 every year. The designation uses the score rounded to the nearest whole point with halves rounded up, so the operative cutoffs on the unrounded score are 89.5 and 79.5. Scores below 60 fail. The same score bands were carried into 24 CFR 5.705(c) under the National Standards for the Physical Inspection of Real Estate (NSPIRE), effective 1 July 2023 for public housing and 1 October 2023 for multifamily programs (88 FR 30442).

The schedule creates two sharp changes in monitoring intensity at arbitrary points of a continuous score. Whether less frequent inspection lets conditions slip is a first-order question for the design of inspection regimes, and one the schedule's own discontinuities can answer.

Two further features motivate the second and third parts of the study. First, the NSPIRE scoring notice (88 FR 43371) assigns a score of 59 to any property whose unit deficiencies account for 30 or more points of deductions, whatever its total, and the published scores show a large mass at exactly 59. Second, HUD's scoring notice of 29 September 2026 (FR Doc. 2026-19881) states that from 1 October 2026 deficiencies under the new affirmative requirements (including fire-labelled doors, GFCI protection, guardrails, HVAC, and interior and minimum electrical and lighting requirements) are included in the score, after two one-year postponements.

### 1.2 Research questions

- **RQ1.** Does assignment to a longer inspection interval change a property's score at its next inspection?
- **RQ2.** Is the excess of NSPIRE scores at exactly 59 fully accounted for by the unit-threshold rule moving properties down from 60 to 70, or is there additional sorting around the passing score? Was there sorting around the 60, 80 and 90 thresholds under the earlier protocol, which had no such rule?
- **RQ3.** Did scores, failure rates and the mass at 59 change when affirmative requirements began to be scored on 1 October 2026?

### 1.3 Hypotheses

**H1 (confirmatory).** For multifamily properties, crossing a schedule cutoff changes the next inspection score.

- H1a: at the 80 cutoff (two-year instead of annual inspection).
- H1b: at the 90 cutoff (three-year instead of two-year inspection).

Tests are two-sided. The directional prediction, if inspection frequency disciplines maintenance, is a lower next score just above each cutoff. A null result is informative: the design can detect about {J['mde_80_h5']:.1f} points at 80 and {J['mde_90_h5']:.1f} points at 90 (Section 6.4).

**H2 (confirmatory).** Under NSPIRE, the excess number of inspections scored exactly 59 equals the number missing from scores 60 to 70, relative to a smooth counterfactual distribution. The unit-threshold rule can only move a property to 59 from a score of 60 to 70, because 30 or more points of unit deductions cap the total at 70. Equality is what a purely mechanical account predicts; an excess at 59 larger than the deficit in 60 to 70, or a deficit elsewhere, indicates sorting beyond the rule.

**H2b (confirmatory, descriptive).** Under the earlier protocol (UPCS), which had no unit-threshold rule, there are more inspections just above each of 60, 80 and 90 than a smooth counterfactual predicts. The technical-review rule in 24 CFR 200.857(d) counts an error as significant when correcting it moves a property across a threshold, which gives owners a reason to contest scores just below one.

**H3 (confirmatory, to be run when data exist).** Among NSPIRE inspections, those conducted on or after 1 October 2026 have (a) lower scores, (b) a higher share below 60 and (c) a higher share at exactly 59 than those conducted before, for the same properties.

## 2. Data

### 2.1 Dataset

Assisted Housing Condition Panel (AHCP) v1.0, doi:10.5281/zenodo.23200319, built by the author from the nine physical-inspection score files HUD publishes. Public, CC BY 4.0. The codebook in the deposit defines every variable. H3 will use a later AHCP release built by the same code once HUD publishes a file containing inspections dated on or after 1 October 2026; AHCP v1.0 ends on {J['latest_date']}.

### 2.2 Samples

**Primary sample (H1).** Multifamily inspections listed in HUD's 2011 history file, which reports up to ten inspections per property for 2000 to 2009 with scores to two decimals. An *index inspection* is one dated from 8 January 2001 (the rule's effective date) to 31 December 2005, which leaves at least 46 months before the file ends in November 2009. There are {n(J['n_index'])} index inspections of {n(J['n_index_props'])} properties; {100 * J['share_has_next']:.1f}% have a later inspection listed. {100 * J['share_noninteger']:.1f}% of index scores are non-integer, so the running variable is effectively continuous.

The index window was fixed before any design check was run and is not changed by them.

**Bunching samples (H2, H2b).** All AHCP inspections flagged NSPIRE, by program ({n(lm['n'])} multifamily, {n(lp['n'])} public housing); and all UPCS inspections dated 2013 to 2019, by program, which carry integer scores.

**Post-change sample (H3).** Properties with at least one NSPIRE inspection on each side of 1 October 2026 in the first AHCP release that contains inspections through at least 31 March 2027.

### 2.3 Why not the later multifamily and public housing data for H1

After 2013 the public files list only each property's latest inspection, so the next inspection is observed for a selected subset, and scores are integers. Public housing adopted the project-level schedule only in 2011, with exceptions for small and troubled agencies that AHCP cannot identify. These samples are used for H2b and for exploratory analysis only.

## 3. Prior knowledge of the data

This is a secondary analysis of data the author assembled. The following is a complete statement of what was known when this plan was written.

**Already published by the author (AHCP data descriptor).** The distribution of scores by protocol, including the mass of NSPIRE scores at exactly 59 and the sparse 60 to 64 band; mean score changes and correlations between consecutive inspections by protocol, pooled over all scores; the distribution of intervals between inspections by period; the direction of score revisions between file vintages.

**Design checks run for this plan (`code/01_design_checks.py`, output in `output/design_checks.txt`).** These use the running variable, the treatment and predetermined variables. They do not compare the outcome across a cutoff.

- Sample sizes near each cutoff. Within 5 points: {near('80', 5)} at 80; {near('90', 5)} at 90.
- Density of the index score at the cutoffs (log ratio of density above to below, local linear on half-point bins, bandwidth 6): {d80['log_ratio']:+.3f} (z = {d80['z']:+.2f}) at 80; {d90['log_ratio']:+.3f} (z = {d90['z']:+.2f}) at 90; {d60['log_ratio']:+.3f} (z = {d60['z']:+.2f}) at 60.
- First stage. Years to the next listed inspection jump by {fs80['jump']:+.2f} (se {fs80['se']:.2f}) at 80, from {fs80['left_limit']:.2f} just below, and by {fs90['jump']:+.2f} (se {fs90['se']:.2f}) at 90, from {fs90['left_limit']:.2f}. The probability of an interval longer than 1.5 years jumps by {pr80['jump']:.2f} at 80; of one longer than 2.5 years, by {pr90['jump']:.2f} at 90.
- Attrition. The probability that a next inspection is listed falls by {abs(at80['jump']):.3f} (se {at80['se']:.3f}) at 80 and {abs(at90['jump']):.3f} (se {at90['se']:.3f}) at 90, and the gap does not shrink when the index window is shortened. {100 * J['attrition_absent_later_share']:.0f}% of index inspections without a next listing belong to properties absent from every later HUD file. The pattern is what a constant exit rate from HUD's portfolio would produce when the next inspection is a year further away. Section 6.3 pre-specifies how it is handled.
- Predetermined covariates at the cutoffs (jump, bandwidth 5): previous score {bal('prev_score_hist', '80')} at 80 and {bal('prev_score_hist', '90')} at 90; inspection sequence number {bal('seq_hist', '80')} and {bal('seq_hist', '90')}; index year {bal('index_year', '80')} and {bal('index_year', '90')}.
- Residual standard deviation of the next score around a quartic polynomial in the index score fitted with no cutoff terms ({J['resid_sd']:.1f} points), used only for the power calculation.
- Counts of integer scores at 59, 60, 79, 80, 89 and 90 in the later samples (for example, {n(lm['59'])} NSPIRE multifamily inspections at 59 and {lm['60']} at 60; {n(um['79'])} UPCS multifamily inspections at 79 and {n(um['80'])} at 80 in 2013 to 2019).

**Not examined.** No contrast of the next score, next failure, or fixed-horizon score across any cutoff has been computed or plotted, at any bandwidth. The excess-mass and missing-mass statistics for H2 and the counterfactual threshold counts for H2b have not been computed on real data. No inspection dated on or after 1 October 2026 exists in the data.

**Safeguards.** The confirmatory scripts (`code/02_rd_analysis.py`, `code/03_bunching_analysis.py`) refuse to run on real outcomes until `docs/REGISTRATION.json` records the public OSF registration. They were tested end to end on a simulated outcome with no discontinuity (`--dry-run`). The estimators were validated on synthetic data with known effects (`code/00_simulation_tests.py`). The repository's commit history timestamps the code that precedes registration.

## 4. Variables

| Role | Variable | Definition |
|---|---|---|
| Running variable | Index score | `inspection_score` of the index inspection, as published (two decimals) |
| Assignment | Above cutoff | Index score at or above 79.5 (H1a) or 89.5 (H1b) |
| Treatment | Interval | Years from the index inspection to the next listed inspection of the same property |
| Primary outcome | Next score | `inspection_score` of the next listed inspection of the same property |
| Secondary outcome | Next failure | Next score, rounded half up, below 60 |
| Secondary outcome | Fixed-horizon score | Score at the first listed inspection at least 21 months after the index inspection |
| Attrition | Has next | A next inspection is listed |
| Predetermined | Previous score, sequence number, index year, index month | From the same history file |
| Cluster | Property | `ahcp_property_id` |

No observations are excluded as outliers. Scores are bounded between 0 and 100. No weights are used.

## 5. Estimands

**Reduced form (primary).** The difference in the expected next score between properties just above and just below a cutoff: the effect of being *assigned* the longer interval. The next score is measured about a year later on the longer-interval side. The estimand therefore answers the policy question as posed, namely what the regulator observes when it next looks, and combines any effect of monitoring with an extra year of elapsed time.

**Fixed-horizon (secondary).** The same contrast for the score at the first inspection at least 21 months after the index inspection. Just below 80 this is typically the second follow-up inspection and just above it the first, so elapsed time is roughly equalized and the contrast is closer to the effect of the intervening inspection.

**Fuzzy (secondary).** The reduced form divided by the first-stage jump in the interval: the effect per additional year between inspections, for properties whose interval is moved by the cutoff.

## 6. Analysis plan for H1

### 6.1 Primary specification

Local linear regression on each side of the cutoff, triangular kernel, bandwidth 5 score points, standard errors clustered by property. The estimate is the coefficient on the above-cutoff indicator in

  y = a + t * 1[x >= c] + b (x - c) + g * 1[x >= c] (x - c) + e,  for |x - c| <= 5,

estimated separately at c = 79.5 and c = 89.5. The bandwidth is fixed in advance so that no tuning choice depends on the outcome; it keeps each window clear of the other cutoff.

### 6.2 Inference

Two-sided tests at the 5% level. H1a and H1b form one family; Holm's correction is applied across the two. Estimates are reported with 95% confidence intervals whatever their significance.

### 6.3 Threats and pre-specified responses

- **Sorting at the cutoff.** The density check does not reject at 90 and is borderline at 80 (z = {d80['z']:+.2f}). The density test is repeated in the confirmatory run, and with `rddensity` if available. Estimates are also reported excluding scores within 0.5 and within 1 point of the cutoff. *Decision rule:* if the density test rejects at the 5% level at a cutoff, the 1-point donut estimate is the headline for that cutoff and is described as resting on extrapolation.
- **Differential attrition.** Because follow-up is less likely above each cutoff, trimming bounds in the manner of Lee (2009) are reported: the below-cutoff outcomes within the bandwidth are trimmed from the top and from the bottom by the excess follow-up share, and the specification is re-estimated. *Decision rule:* a primary result is described as robust only if the bounds exclude zero; the 95% interval for the bounds is reported.
- **Covariate imbalance.** The four predetermined variables are tested at each cutoff with the primary specification.
- **Specification.** Bandwidths of 3 and 8; a local quadratic at bandwidth 8; placebo cutoffs at 74.5, 84.5 and 94.5 with bandwidth 4; exclusion of each property's first listed inspection, whose own schedule position is unknown.
- **Data-driven bandwidth.** The primary estimates are also computed with an MSE-optimal bandwidth and robust bias-corrected confidence intervals (Calonico, Cattaneo and Titiunik 2014) using `rdrobust`, clustered by property. If the software cannot be installed this is reported as a deviation.

### 6.4 Power

With the primary specification and the residual standard deviation of {J['resid_sd']:.1f} points, the simulated standard error of the reduced-form jump is {J['mde_80_h5'] / 2.8:.2f} at 80 and {J['mde_90_h5'] / 2.8:.2f} at 90. The minimum detectable effects at 80% power and a 5% two-sided test are {J['mde_80_h5']:.1f} and {J['mde_90_h5']:.1f} points, about {J['mde_80_h5'] / J['resid_sd']:.2f} and {J['mde_90_h5'] / J['resid_sd']:.2f} residual standard deviations. Effects smaller than this cannot be ruled out by a null result, and the confidence interval will be reported as the range of effects the data support.

## 7. Analysis plan for H2 and H2b

**H2.** For NSPIRE inspections, separately for multifamily and public housing, the counts at each integer score from 40 to 85 are fitted with a fifth-order polynomial in the score (log counts), excluding scores 55 to 72. The excess at 59 is the observed count minus the fitted count; the missing mass is the fitted minus observed count summed over 60 to 70. The test statistic is excess minus missing, with a standard error from 500 Poisson resamples of the counts. Holm's correction is applied across the two programs. Robustness: polynomial orders 4 and 6; excluded windows 53 to 74 and 56 to 71. The same statistics are reported for UPCS inspections from 2013 to 2019 as a comparison in which no excess at 59 is expected.

*Interpretation rule.* If excess minus missing is not distinguishable from zero and the excess is large, the mass at 59 is described as consistent with the unit-threshold rule alone. The public files do not report unit deductions, so the rule's application cannot be observed directly; the conclusion is about consistency with the mechanical account, not proof of it.

**H2b.** For UPCS inspections from 2013 to 2019, by program, and for each threshold c in 60, 80 and 90: counts at integer scores from c - 10 to c + 9 are fitted with a quadratic excluding c - 2 to c + 1, and observed and fitted counts in the two scores below and the two at or above the threshold are compared. For the 2001 to 2005 history sample, the density discontinuity at 59.5 is reported with the same estimator as the design check (already seen: z = {d60['z']:+.2f}).

## 8. Analysis plan for H3

**Trigger.** The analysis is run once, on the first HUD score file whose inspections extend to at least 31 March 2027, added to AHCP by the released pipeline.

**Specification.** Among NSPIRE inspections of properties with at least one NSPIRE inspection on each side of 1 October 2026, regress each outcome (score; score below 60; score exactly 59) on an indicator for inspection on or after 1 October 2026, property fixed effects, and program-by-inspection-year effects, clustering by property. Holm's correction across the three outcomes.

**Limits stated in advance.** The public files list only each property's latest inspection per vintage, so a property contributes a before-and-after pair only if an earlier vintage captured the earlier inspection. Properties reinspected soon after the change are disproportionately those with low earlier scores, which biases a naive comparison toward improvement; the property fixed effects do not remove this. The result will be reported as an association with that caveat, alongside the unadjusted distributions before and after.

## 9. Exploratory analyses

Labelled exploratory in any paper: the effect of failing at the 60 cutoff on the next score, where sorting is already evident; heterogeneity of H1 by previous score and by region; H1 in the 2013 to 2019 multifamily sample with a discrete running variable; public housing after 2011; persistence of effects at the second next inspection.

## 10. Deviations and reporting

Any departure from this plan will be listed in a table of deviations with its reason. All pre-specified results are reported, including nulls. Code, the design-check output and the results files are public in the project repository, and the confirmatory outputs are produced by the scripts named above without manual steps.

## 11. Software

Python 3 with NumPy, pandas and SciPy; estimators in `code/rdlib.py`; `rdrobust` and `rddensity` for the robustness analyses in Section 6.3.

## References

- Calonico, S., Cattaneo, M. D., and Titiunik, R. (2014). Robust nonparametric confidence intervals for regression-discontinuity designs. *Econometrica*, 82(6), 2295-2326.
- Lee, D. S. (2009). Training, wages, and sample selection: estimating sharp bounds on treatment effects. *Review of Economic Studies*, 76(3), 1071-1102.
- Lee, D. S., and Lemieux, T. (2010). Regression discontinuity designs in economics. *Journal of Economic Literature*, 48(2), 281-355.
- McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: a density test. *Journal of Econometrics*, 142(2), 698-714.
- Kleven, H. J. (2016). Bunching. *Annual Review of Economics*, 8, 435-464.
- Ameen, A. (2026). Assisted Housing Condition Panel (AHCP) v1.0. Zenodo. doi:10.5281/zenodo.23200319.
- U.S. Department of Housing and Urban Development. Uniform physical condition standards and physical inspection requirements for certain HUD housing; administrative process for assessment of insured and assisted properties. Final rule, 65 FR 77230, 8 December 2000.
- U.S. Department of Housing and Urban Development. Economic Growth Regulatory Relief and Consumer Protection Act: implementation of NSPIRE. Final rule, 88 FR 30442, 11 May 2023.
- U.S. Department of Housing and Urban Development. NSPIRE and associated protocols, scoring notice. 88 FR 43371, 7 July 2023.
- U.S. Department of Housing and Urban Development. National Standards for the Physical Inspection of Real Estate: scoring notice; request for public comment. FR Doc. 2026-19881, 29 September 2026.
"""
os.makedirs(os.path.join(ROOT, "docs"), exist_ok=True)
open(os.path.join(ROOT, "docs", "PAP.md"), "w", encoding="utf-8").write(doc); print("PAP written", len(doc.split()), "words")
