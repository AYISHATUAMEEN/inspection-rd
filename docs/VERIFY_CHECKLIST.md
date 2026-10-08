# Verification checklist

## Before registering (now)
- [ ] `docs/PAP.md` read in full; hypotheses, primary specification (local linear, bandwidth 5, clustered by property) and decision rules are ones you will stand behind
- [ ] Section 3 is a complete and true statement of what you have seen of the data
- [ ] Regulatory facts in section 1.1 checked against the cited rules: 24 CFR 200.857(b) score bands and rounding; effective date 8 January 2001; NSPIRE dates; the 59 unit-threshold rule; the 1 October 2026 scoring change
- [ ] `python code/00_simulation_tests.py` passes on your machine
- [ ] `python code/01_design_checks.py` reproduces `output/design_checks.txt`
- [ ] `python code/02_rd_analysis.py` without `--dry-run` is refused (the gate works)

## After the confirmatory run
- [ ] Every number in the paper traced to `output/rd_results.json` or `output/bunching_results.json`
- [ ] Deviations table complete
- [ ] Paper read and revised in your own voice
