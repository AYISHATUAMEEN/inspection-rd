# Do Inspection Thresholds Improve Housing Conditions?

**Status:** pre-registration stage. The pre-analysis plan is written; no outcome has been analysed. Confirmatory scripts are locked until the plan is registered on OSF.

HUD inspects assisted multifamily properties every three years if they score 90 or above, every two years at 80 to 89, and every year below 80. This project uses those cutoffs in a regression discontinuity design to estimate whether a longer inspection interval changes a property's next score, tests whether the mass of NSPIRE scores at exactly 59 is accounted for by the unit-threshold scoring rule, and pre-specifies an analysis of the 1 October 2026 change that began scoring affirmative requirements.

Data: Assisted Housing Condition Panel v1.0, doi:10.5281/zenodo.23200319. Place `ahcp_inspections.csv.gz` and `ahcp_properties.csv` in `data/input/`.

## What exists now

| File | Content |
|---|---|
| `docs/PAP.md`, `docs/PAP.pdf` | Pre-analysis plan |
| `output/design_checks.txt` | Sample sizes, density, first stage, attrition, balance, power. No outcome contrasts |
| `code/rdlib.py` | Local-linear RD (sharp and fuzzy), density ratio, bunching estimator |
| `code/00_simulation_tests.py` | Validation of the estimators on synthetic data |
| `code/02_rd_analysis.py`, `code/03_bunching_analysis.py` | Confirmatory analyses; refuse to run before registration; `--dry-run` uses a simulated outcome |
| `docs/OSF_GUIDE.md` | How to register the plan |

## Design in brief

54,915 index inspections of 28,936 multifamily properties, 2001 to 2005. At the 80 cutoff the interval to the next inspection jumps by 1.08 years; at 90 by 0.99 years. Minimum detectable effects on the next score are about 1.5 and 1.4 points.

## Run

```
python code/00_simulation_tests.py
python code/01_design_checks.py
python code/04_write_plan.py
python code/02_rd_analysis.py --dry-run      # real run only after registration
python code/03_bunching_analysis.py --dry-run
```

Code: MIT. Documents: CC BY 4.0. Maintainer: Ayishatu Ameen (ameenayishatu4@gmail.com).

AI assistance: code and drafts were prepared with Claude (Anthropic); the author makes the analytic decisions and verifies the outputs.
