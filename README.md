# Associations of vaccination and monoclonal antibody therapy with long COVID and mortality

Statistical analysis code accompanying the manuscript examining the effects of pre-infection
mRNA vaccination status and monoclonal antibody treatment on post-acute sequelae of
SARS-CoV-2 infection (PASC), disease severity, and 3-year all-cause mortality.


## Overview

The analysis estimates average treatment effects (ATE) in an entropy-balanced
pseudo-population. Two exposures are evaluated using independent weighting models:

- **mRNA vaccination status** (multi-category): unvaccinated, fully vaccinated/boosted over
  6 months, fully vaccinated/boosted within 6 months
- **Monoclonal antibody treatment** (binary)

Covariate selection for all propensity models follows the directed acyclic graph described in
the manuscript Methods:

| Model | Covariates |
|---|---|
| Vaccination | age, sex, race, ethnicity, Charlson Comorbidity Index, prior infection |
| Monoclonal antibody | the above, plus disease severity and vaccination status |

Post-infection variables (disease severity, monoclonal antibody, antiviral therapy) are
excluded from vaccination models as they lie on the causal pathway from vaccination to PASC.
Antiviral therapy is excluded from monoclonal antibody models because both treatments share
similar indications and adjustment could induce collider bias.

## Repository contents

| File | Description |
|---|---|
| `Causal_Analysis_Primary.R` | Entropy-balanced weighted logistic models for the effect of vaccination (`outglm`) and monoclonal antibody (`outglm_mAB`) on PASC. Used for the primary analysis, the prespecified age-stratified analyses (18-40, 40-80, >80), and disease severity as a secondary outcome. |
| `Survival_Analysis.R` | Entropy-balanced weighted Cox proportional hazards models for 3-year all-cause mortality, with proportional hazards assessed via scaled Schoenfeld residuals. |
| `Sensitivity_Analysis.R` | Active comparator analysis restricting to patients receiving monoclonal antibody or antiviral therapy exclusively (`sens_mab_confounding`), and antiviral therapy as a negative control exposure (`sens_antiviral_control`). |

## Scope

**These scripts contain model specifications only.** Cohort construction and data processing
are not included, as they depend on institutional data structures that cannot be shared.
Each script expects an analytic dataset (`dat`) that has already been assembled with the
following:

- One row per patient, restricted to adults with confirmed SARS-CoV-2 infection
- Partially vaccinated individuals (single mRNA dose prior to infection) excluded
- `vaccination_status_cat` (factor: 0 = unvaccinated, 2 = vaccinated/boosted over 6 months,
  3 = vaccinated/boosted within 6 months)
- `received_mAb`, `anti.viral`, `prior_infection` (binary)
- `age`, `sex_cd`, `race`, `hispanic`, `CHARLSON_INDEX`, `severity`, `severity_binary`
- `PASC.any` (outcome)
- For the survival script: `LabStartDTS` (index infection date) and `death_date`, with
  patients who died within 28 days of infection excluded

## Requirements

R version 4.3.0 or later.

```r
install.packages(c("dplyr", "tidyr", "readr", "WeightIt", "cobalt", "marginaleffects",
                   "survival", "survey", "survminer", "broom", "tidycmprsk", "ggsurvfit",
                   "ggplot2", "flextable", "officer", "gtsummary"))
```

Entropy balancing weights are estimated with **WeightIt**, balance is assessed with
**cobalt**, marginal contrasts are computed with **marginaleffects**, and weighted survival
models use **survival** and **survey**.

## Usage

Replace the `dat` placeholder at the top of each script with your analytic dataset, then run.

```r
# Primary analysis: vaccination effect on PASC
outglm("PASC.any", dat, group = "all", severity = F)

# Primary analysis: monoclonal antibody effect on PASC
outglm_mAB("PASC.any", dat, group = "all", vax = T, severity = T)

# Active comparator sensitivity analysis
sens_mab_confounding(dat.aoi = dat, group = "Overall")
```

Balance should be checked before interpreting any effect estimate. The `bal.tab()` calls in
each script report standardized mean differences, variance ratios, and KS statistics.

## License

MIT — see `LICENSE` for details.