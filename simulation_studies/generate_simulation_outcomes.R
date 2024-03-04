###########################################
## function for simulations
## data generating model
## generate the outcomes
## this function takes a list of object as input, which was produced by the function "simulate_component_data()"
###########################################


simulate_outcomes_few_scenarios <- function(dat, True_Beta = 1.5, ancestry_var_eff = 0.5){
  
  ###########################################
  ## homogeneous weighting PRS setting
  ###########################################
  ## PRS + conf-PRS1
  O1 <- generate_outcome(dat$standard_prs, 
                         ancestry_var = dat$conf_PRS1, 
                         True_Beta,
                         ancestry_var_eff)
  
  ## PRS + conf-PRS2
  O2 <- generate_outcome(dat$standard_prs, 
                         ancestry_var = dat$conf_PRS2, 
                         True_Beta,
                         ancestry_var_eff)
  
  ## PRS + single variant
  O3 <- generate_outcome(dat$standard_prs, 
                         ancestry_var = dat$single_conf_var, 
                         True_Beta,
                         ancestry_var_eff)
  
  ## PRS + two variants
  ## change to generate_outcome2() function
  ## because the model includes two variables as covariates
  O4 <- generate_outcome2(dat$standard_prs, 
                         ancestry_var = dat$two_conf_vars, 
                         True_Beta,
                         ancestry_var_eff)
  ## PRS + conf-pc*
  O5 <- generate_outcome(dat$standard_prs, 
                         ancestry_var = dat$conf_pc, 
                         True_Beta,
                         ancestry_var_eff)
  
  ## PRS (no-confounding)
  O6 <- generate_outcome(dat$standard_prs, 
                         ancestry_var = NULL, 
                         True_Beta,
                         ancestry_var_eff)
  
  ###########################################
  ## heterogeneous weighting PRS setting
  ###########################################
  ## PRS + conf-PRS1
  O7 <- generate_outcome(dat$ancestry_effect_prs, 
                         ancestry_var = dat$conf_PRS1, 
                         True_Beta,
                         ancestry_var_eff)
  
  ## PRS + conf-PRS2
  O8 <- generate_outcome(dat$ancestry_effect_prs, 
                         ancestry_var = dat$conf_PRS2, 
                         True_Beta,
                         ancestry_var_eff)
  
  ## PRS + single variant
  O9 <- generate_outcome(dat$ancestry_effect_prs, 
                         ancestry_var = dat$single_conf_var, 
                         True_Beta,
                         ancestry_var_eff)
  
  ## PRS + two variants
  O10 <- generate_outcome2(dat$ancestry_effect_prs, 
                           ancestry_var = dat$two_conf_vars,
                           True_Beta,
                           ancestry_var_eff)
  
  ## PRS + conf-pc*
  O11 <- generate_outcome(dat$ancestry_effect_prs,
                          ancestry_var = dat$conf_pc,
                          True_Beta,
                          ancestry_var_eff)
  
  ## PRS (no-confounding)
  O12 <- generate_outcome(dat$ancestry_effect_prs,
                          ancestry_var = NULL,
                          True_Beta,
                          ancestry_var_eff)
  
  
  return(list(ancestry_var_eff = ancestry_var_eff,
              mod_sPRS_conf_PRS1 = O1,
              mod_sPRS_conf_PRS2 = O2,
              mod_sPRS_conf_singleVar = O3,
              mod_sPRS_conf_twoVars = O4,
              mod_sPRS_conf_pc = O5,
              mod_sPRS_conf_none = O6,
              mod_ancPRS_conf_PRS1 = O7,
              mod_ancPRS_conf_PRS2 = O8,
              mod_ancPRS_conf_singleVar = O9,
              mod_ancPRS_conf_twoVars = O10,
              mod_ancPRS_conf_pc = O11,
              mod_ancPRS_conf_none = O12))
}




generate_outcome <- function(prs, ancestry_var, True_Beta, ancestry_var_eff){
  # the effect is in SD of the ancestry_var, so standardize ancestry_var
  if (is.null(ancestry_var)){
    outcome <- 1 + (True_Beta*prs) + rnorm(length(prs))
    return(outcome)
  } 
  outcome <- 1 + (True_Beta*prs) + 
                (ancestry_var/sd(ancestry_var))*ancestry_var_eff + 
                rnorm(length(prs))
  return(outcome)
    
} 


## for two-variants confounding setting
generate_outcome2 <- function(prs, ancestry_var, True_Beta, ancestry_var_eff){
  outcome <- 1 + (True_Beta*prs) + 
    (ancestry_var[,1]/sd(ancestry_var[,1]))*ancestry_var_eff + 
    (ancestry_var[,2]/sd(ancestry_var[,2]))*ancestry_var_eff + rnorm(length(prs))
    return(outcome)
}


