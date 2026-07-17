library(WeightIt)
library(marginaleffects)
library(cobalt)
library(dplyr)

# Import your data
dat <- "Your study data"

# --- ACTIVE COMPARATOR: mAb vs Antiviral ---
# Restricted to patients receiving either treatment exclusively

sens_mab_confounding <- function(dat.aoi, group, severity = T) {
  
  # Restrict to patients who received EITHER mAb OR antiviral (exclusive)
  # This compares two treatments given to similar high-risk populations
  dat.aoi <- dat.aoi %>%
    filter((received_mAb == 1 & anti.viral == 0) | (received_mAb == 0 & anti.viral == 1)) %>%
    mutate(
      treatment = factor(ifelse(received_mAb == 1, "mAb", "Antiviral"), 
                         levels = c("Antiviral", "mAb"))
    )
  dat.aoi$label <- as.factor(dat.aoi$PASC.any)
  dat.aoi$label <- ifelse(dat.aoi$label == "0",0,1)#as.numeric(dat.aoi$label)
  
  cat("\n=== Sample Sizes ===\n")
  cat("Antiviral only:", sum(dat.aoi$treatment == "Antiviral"), "\n")
  cat("mAb only:", sum(dat.aoi$treatment == "mAb"), "\n")
  cat("Long COVID events:", sum(dat.aoi$label), "\n\n")
  
  # Entropy balancing
  if (severity == T) {
    W.out <- weightit(
      treatment ~ age + sex_cd + CHARLSON_INDEX + severity + race + 
        hispanic  + vaccination_status_cat +prior_infection,
      data = dat.aoi,
      estimand = "ATE",
      method = "ebal"
    )
  } else {
    W.out <- weightit(
      treatment ~ age + sex_cd + CHARLSON_INDEX  + race + 
        hispanic  + vaccination_status_cat +prior_infection ,
      data = dat.aoi,
      estimand = "ATE",
      method = "ebal"
    )
  }
  
  cat("=== Balance ===\n")
  bal <- bal.tab(W.out, stats = c("m", "v", "ks"), m.threshold = .05, disp.v.ratio = TRUE)
  print(bal)
  
  # Weighted logistic regression
  fit <- glm_weightit(
    label ~ treatment,
    data = dat.aoi,
    weightit = W.out,
    family = binomial(logit)
  )
  
  # Odds ratio
  output_or <- avg_comparisons(fit, variables = "treatment", 
                               comparison = "lnoravg", transform = "exp")
  
  # Risk difference
  output_rd <- avg_comparisons(fit, variables = "treatment", 
                               comparison = "difference")
  
  results <- data.frame(
    analysis = "mAb_vs_Antiviral",
    group = group,
    n_mab = sum(dat.aoi$treatment == "mAb"),
    n_antiviral = sum(dat.aoi$treatment == "Antiviral"),
    OR = output_or$estimate,
    OR_lower = output_or$conf.low,
    OR_upper = output_or$conf.high,
    OR_p = output_or$p.value,
    RD = output_rd$estimate,
    RD_lower = output_rd$conf.low,
    RD_upper = output_rd$conf.high
  )
  
  return(results)
}


# --- SUPPORTING ANALYSIS: Antiviral as Negative Control Exposure ---
# If antiviral shows similar "harmful" effect, confirms confounding by indication

sens_antiviral_control <- function(dat.aoi, group, severity = T) {
  
  dat.aoi$label <- as.factor(dat.aoi$PASC.any)
  dat.aoi$label <- ifelse(dat.aoi$label == "0",0,1)#as.numeric(dat.aoi$label)
  
  cat("\n=== Antiviral Distribution ===\n")
  cat("No antiviral:", sum(dat.aoi$anti.viral == 0), "\n")
  cat("Antiviral:", sum(dat.aoi$anti.viral == 1), "\n\n")
  
  if (severity == T) {
    W.out <- weightit(
      anti.viral ~ age + sex_cd + CHARLSON_INDEX + severity + prior_infection + race + 
        hispanic  + vaccination_status_cat,
      data = dat.aoi,
      estimand = "ATE",
      method = "ebal"
    )
  } else {
    W.out <- weightit(
      anti.viral ~ age + sex_cd + CHARLSON_INDEX + prior_infection + race + 
        hispanic  + vaccination_status_cat,
      data = dat.aoi,
      estimand = "ATE",
      method = "ebal"
    )
  }
  
  cat("=== Balance ===\n")
  bal <- bal.tab(W.out, stats = c("m", "v", "ks"), m.threshold = .05, disp.v.ratio = TRUE)
  print(bal)
  
  fit <- glm_weightit(
    label ~ anti.viral,
    data = dat.aoi,
    weightit = W.out,
    family = binomial(logit)
  )
  
  output_or <- avg_comparisons(fit, variables = "anti.viral", 
                               comparison = "lnoravg", transform = "exp")
  output_rd <- avg_comparisons(fit, variables = "anti.viral", 
                               comparison = "difference")
  
  results <- data.frame(
    analysis = "Antiviral_vs_None",
    group = group,
    n_treated = sum(dat.aoi$anti.viral == 1),
    n_untreated = sum(dat.aoi$anti.viral == 0),
    OR = output_or$estimate,
    OR_lower = output_or$conf.low,
    OR_upper = output_or$conf.high,
    OR_p = output_or$p.value,
    RD = output_rd$estimate,
    RD_lower = output_rd$conf.low,
    RD_upper = output_rd$conf.high
  )
  
  return(results)
}

# Examples
# Main sensitivity analysis
results_active_comparator <- sens_mab_confounding(dat.aoi = dat, group = "Overall")
results_active_comparator

# Negative control analysis
results_antiviral <- sens_antiviral_control(dat.aoi = dat, group = "Overall")
results_antiviral