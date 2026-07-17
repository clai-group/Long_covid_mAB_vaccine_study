pacman::p_load(tidyverse, survival, survey, survminer, WeightIt, cobalt, 
               gtsummary, broom, tidycmprsk, ggsurvfit)

study_df_mort = "your study dat"

# PS WEIGHTING
# VACCINE
ps_vax_mort <- weightit(
  vaccination_status_cat ~ age + sex_cd + race + hispanic + CHARLSON_INDEX + 
    prior_infection ,
  data = study_df_mort,
  method = "ebal",
  estimand = "ATE"
)

study_df_mort$ps_weights_vax <- ps_vax_mort$weights
bal.tab(ps_vax_mort, stats = c("m", "v"), m.threshold = 0.1)

# mAB
ps_mab_mort <- weightit(
  received_mAb ~ age + sex_cd + race + hispanic + CHARLSON_INDEX + 
    severity + prior_infection + vaccination_status_cat,
  data = study_df_mort,
  method = "ebal",
  estimand = "ATE"
)

study_df_mort$ps_weights_mab <- ps_mab_mort$weights
bal.tab(ps_mab_mort, stats = c("m", "v"), m.threshold = 0.1)



# COX  
# --- Vaccination ---
design_vax_mort <- svydesign(ids = ~1, weights = ~ps_weights_vax, data = study_df_mort)

cox_vax_mort <- svycoxph(
  Surv(follow_up_days, death_3yr) ~ vaccination_status_cat,
  design = design_vax_mort
)

cox_vac_tbl = tidy(cox_vax_mort, conf.int = TRUE, exponentiate = TRUE)
cox_vac_tbl
# --- Monoclonal Antibody ---
design_mab_mort <- svydesign(ids = ~1, weights = ~ps_weights_mab, data = study_df_mort)

cox_mab_mort <- svycoxph(
  Surv(follow_up_days, death_3yr) ~ received_mAb,
  design = design_mab_mort
)

cox_mab_tbl = tidy(cox_mab_mort, conf.int = TRUE, exponentiate = TRUE)
cox_mab_tbl

# PH ASSUMPTION
cox_vax_unweighted <- coxph(Surv(follow_up_days, death_3yr) ~ vaccination_status_cat, 
                            data = study_df_mort)
cox.zph(cox_vax_unweighted)

cox_mab_unweighted <- coxph(Surv(follow_up_days, death_3yr) ~ received_mAb, 
                            data = study_df_mort)
cox.zph(cox_mab_unweighted)

