# Registering the plan on OSF (about 20 minutes)

Registration is your attestation that the plan precedes the analysis, so you submit it yourself. A registration cannot be edited or deleted once submitted; it can only be withdrawn, and the withdrawal stays visible.

## Before you start

- Read `docs/PAP.md` end to end. Change anything you disagree with now; tell Claude so the code is changed to match before you register.
- Check section 3 ("Prior knowledge of the data"). It must be a complete and true account of what you have seen.

## Steps

1. Sign in at https://osf.io (ORCID sign-in works). Create a project: **Create new project**, title *Do Inspection Thresholds Improve Housing Conditions?*, storage location of your choice.
2. In the project, **Files**: upload `docs/PAP.pdf`, `output/design_checks.txt`, and a zip of the `code/` folder. These are the materials the registration freezes.
3. **Registrations** tab, **New registration**, template **Secondary Data Preregistration**.
4. Fill the form from the plan:

   | OSF form section | Paste from `docs/PAP.md` |
   |---|---|
   | Title, contributors | Title; you as sole author |
   | Study information: research questions, hypotheses | Sections 1.2 and 1.3 |
   | Data description: datasets, availability, access, date, collection | Section 2.1 and 2.2; public, doi:10.5281/zenodo.23200319 |
   | Variables: manipulated, measured, units, missing data, outliers, weights | Sections 4 and 5 |
   | Knowledge of data: prior publications, prior knowledge | Section 3, in full |
   | Analyses: models, effect size, power, inference criteria, assumption violations, robustness, exploratory | Sections 6 to 9 |
   | Other | Section 10; attach `PAP.pdf` as the authoritative text |

5. Choose **Make registration public immediately** (an embargo of up to four years is available if you prefer; either satisfies the gate).
6. Submit and approve the confirmation email.
7. Copy the registration URL (it looks like `https://osf.io/abcde`) and its DOI. Send both to Claude, or create `docs/REGISTRATION.json` yourself from `docs/REGISTRATION.example.json`.

After that the confirmatory scripts will run.
