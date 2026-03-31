library(dplyr)
library(readr)
library(cobalt)
library(marginaleffects)
library(WeightIt)
library(flextable)
library(ggplot2)
library(officer)
#library(mediation)
library(lme4)
library(DescTools)
library(tidyr)

# Import your data
dat <- "Your study data"
# Covariates selected based on DAG (pre-infection confounders)
covariates_vax <- c("age", "sex_cd", "race", "hispanic", "CHARLSON_INDEX")

# Additional variables for mAb analysis
covariates_mab <- c(
  covariates_vax, "severity", "vaccination_status_cat"
)

### VAcccination
outglm <- function(outcome,#outcome of interest
                   dat.aoi,#data for modeling
                   group,
                   vax = F, # stratify based on this?
                   severity = T,
                   exposure = "vaccination_status_cat",
                   all = T
                   
){
  dat.aoi$label <- as.factor(dat.aoi[[outcome]])
  dat.aoi$label <- ifelse(dat.aoi$label == "0",0,1)#as.numeric(dat.aoi$label)
  
  if(all == T){
    if(severity == T){
      
      dat.aoi <- dat.aoi %>% dplyr::select(all_of(c("vaccination_status_cat","severity", covariates_vax, "label" ))) %>%
        na.omit()
      W.out <- weightit(vaccination_status_cat ~age+sex_cd+CHARLSON_INDEX + severity +
                          race +hispanic +anti.viral,
                        data = dat.aoi, estimand = "ATE" , method = "ebal")
    } else {
      dat.aoi <- dat.aoi %>% dplyr::select(all_of(c("vaccination_status_cat", covariates_vax, "label" ))) %>%
        na.omit()
      print("hi")
      W.out <- weightit(vaccination_status_cat ~ age+sex_cd+CHARLSON_INDEX +race +
                          hispanic +anti.viral,
                        data = dat.aoi, estimand = "ATE" , method = "ebal")
    }
  } 
  bal <- bal.tab(W.out, stats = c("m", "v", "ks"), m.threshold = .05, disp.v.ratio = TRUE,poly = 1)
  #print(bal)
  fit = glm_weightit(label ~ vaccination_status_cat,
                     data = dat.aoi, weightit = W.out,family = binomial(logit) )
  output = avg_comparisons(fit,
                           variables = "vaccination_status_cat",
                           comparison = "lnoravg",
                           transform = "exp")
  
  output = data.frame(output, group)
  
  return(output)
  
  
}

# Example
outglm("PASC.any",dat,group="all")


### mAB
outglm_mAB <- function(outcome,#outcome of interest
                       dat.aoi,#data for modeling
                       group,
                       vax = F, # stratify based on this?
                       severity = T,
                       exposure = "received_mAb",
                       all = T
                       
){
  dat.aoi$label <- as.factor(dat.aoi[[outcome]])
  dat.aoi$label <- ifelse(dat.aoi$label == "0",0,1)#as.numeric(dat.aoi$label)
  
  if (all == TRUE) {
    
    if (severity == TRUE) {
      
      dat.aoi <- dat.aoi %>% dplyr::select(all_of(c("received_mAb",covariates_mab, "label" ))) %>%
        na.omit()
      
      covars <- "age + sex_cd + CHARLSON_INDEX + severity + race + hispanic  +anti.viral"
      
      if (vax == TRUE) {
        covars <- paste(covars, "+ vaccination_status_cat")
      }
      
      fml <- as.formula(paste("received_mAb ~", covars))
      
      w.out <- weightit(
        formula   = fml,
        data      = dat.aoi,
        estimand  = "ATE",
        method    = "ebal"
      )
      
    } else {
      
      dat.aoi <- dat.aoi %>% dplyr::select(all_of(c("received_mAb",covariates_mab, "label" ))) %>%
        na.omit()
      
      covars <- "age + sex_cd + CHARLSON_INDEX + race + hispanic+anti.viral"
      
      if (vax == TRUE) {
        covars <- paste(covars, "+ vaccination_status_cat")
      }
      
      fml <- as.formula(paste("received_mAb ~", covars))
      
      w.out <- weightit(
        formula   = fml,
        data      = dat.aoi,
        estimand  = "ATE",
        method    = "ebal"
      )
    }
  }
  
  bal <- bal.tab(w.out, stats = c("m", "v", "ks"), m.threshold = .05, disp.v.ratio = TRUE,poly = 1)
  #print(bal)
  fit = glm_weightit(label ~ received_mAb,
                     data = dat.aoi, weightit = w.out,family = binomial(logit) )
  output = avg_comparisons(fit,
                           variables = "received_mAb",
                           comparison = "lnoravg",
                           transform = "exp")
  
  output = data.frame(output, group)
  
  return(output)
  
  
}

# Example
outglm_mAB("PASC.any",dat,group="all", vax = T, severity = T)
