###########################################
## function for simulations
## function that estimates PRS associations across the different scenarios of unknown confounding and of PRS
## based on different association models
###########################################

estimate_prs_associations <- function(prs, 
                                      outcomes,
                                      conf_pc,
                                      Data.PC10,
                                      Data.PC20,
                                      proportions_ancestry,
                                      gePRS,
                                      lePRS){
  
  true_confounder <- c(PRS_Conf1 = "PRS1",
                       PRS_Conf2 = "PRS2", 
                       singelVar = "singleVar",
                       twoVars = "twoVars",
                       conf_pc = "pc",
                       none = "none")
  
  true_prs <- c(homogeneous_PRS = "sPRS", 
                heterogeneous_PRS = "ancPRS")
  
  analysis_adjustment <- c("none", "conf_pc", "Data.PC10", "Data.PC20", "gePRS", "lePRS", "proportion_ancestry")
 
  res <- c()
  for (i in 1:length(true_confounder)){
    for (j in 1:length(true_prs)){
      for (k in 1:length(analysis_adjustment)){

        cur_out <- outcomes[[paste0("mod_", true_prs[[j]], "_conf_",  true_confounder[[i]])]]
        
        if (analysis_adjustment[k] == "none"){
          cur_mod <- lm(cur_out ~ prs)  
        }
        if (analysis_adjustment[k] == "conf_pc"){
          cur_mod <- lm(cur_out ~ prs + conf_pc)  
        }
        
        if (analysis_adjustment[k] == "Data.PC10"){
          cur_mod <- lm(cur_out ~ prs + Data.PC10)  
        }
        
        if (analysis_adjustment[k] == "Data.PC20"){
          cur_mod <- lm(cur_out ~ prs + Data.PC20)  
        }
        
        if (analysis_adjustment[k] == "gePRS"){
          grPRS <- prs - gePRS
          cur_mod <- lm(cur_out ~ grPRS + gePRS)  
        }
        if (analysis_adjustment[k] == "lePRS"){
          lrPRS <- prs - lePRS
          cur_mod <- lm(cur_out ~ lrPRS + lePRS)  
        }
        if (analysis_adjustment[k] == "proportion_ancestry"){
          cur_mod <- lm(cur_out ~ prs + as.matrix(proportions_ancestry))
        }
        
        est_coef <- summary(cur_mod)$coef
        cur_res <- data.frame(true_confounder = names(true_confounder)[i],
                              true_prs = names(true_prs)[j],
                              model_adjustment = analysis_adjustment[k],
                              beta = est_coef[2,"Estimate"], 
                              SE = est_coef[2,"Std. Error"], 
                              pval = est_coef[2,"Pr(>|t|)"])

        
        res <- rbind(res, cur_res)
        
      }
    }
  }
  return(res)
  
}
