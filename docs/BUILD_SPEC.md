# Build spec

```
PROJECT:        Do Inspection Thresholds Improve Housing Conditions? Archetypes: pre-registered analysis paper
                (regression discontinuity) with a replication package.
QUESTION:       Does being inspected less often make a HUD-assisted property score worse at its next inspection,
                and is the pile-up of scores at exactly 59 produced by a scoring rule or by behaviour?
SOURCES:        Assisted Housing Condition Panel v1.0, doi:10.5281/zenodo.23200319 (CC BY 4.0), built from HUD's
                public physical-inspection score files. Regulations: 24 CFR 200.857 (65 FR 77230), 24 CFR 5.705,
                88 FR 30442, 88 FR 43371, FR Doc. 2026-19881.
UNIT:           Inspection (index inspection paired with the property's next inspection); property is the cluster.
MEASURES:       Running variable: index score. Cutoffs 79.5 and 89.5. Treatment: years to next inspection.
                Outcome: next score. Bunching: excess at 59 and missing mass in 60-70 against a polynomial
                counterfactual. Full definitions in docs/PAP.md.
OUTPUTS:        docs/PAP.md and PAP.pdf (pre-analysis plan); output/design_checks.txt; after registration:
                output/rd_results.*, output/bunching_results.*, the paper, and SSRN and MPRA submission kits.
VENUES:         OSF Registries (plan, first) -> GitHub + Zenodo (replication package) -> SSRN and MPRA (paper).
VERIFY POINTS:  (1) the plan's hypotheses, primary specification and decision rules; (2) the disclosure of prior
                knowledge in PAP section 3; (3) the regulatory facts in PAP section 1.1; (4) after analysis, every
                reported number against the output files.
LICENSE:        Code MIT; documents CC BY 4.0.
ASSUMPTIONS:    Author block reused from the AHCP project in the same session. The 2011 history file is treated as
                a complete listing of 2000-2009 inspections. H3 cannot be run until HUD publishes inspections dated
                after 1 October 2026. No R package is planned: the replication package is Python.
```

## Order of work

1. Design checks and plan (done before registration).
2. **Gate: the author registers the plan on OSF.** Nothing below runs before this.
3. Confirmatory analysis, figures, paper.
4. Author verification, then SSRN and MPRA submission by the author.
5. H3 when the data exist.
