
# generate EPRS for the standard PRS: weighted sum of alleles 
# (further divided by the number of SNPs used), 
# according to local genetic ancestry

generate_local_EPRS_rPRS <- function(prs, 
                                      eafs, 
                                      betas_vec,
                                      local_anc_counts){
  stopifnot(length(betas_vec) == nrow(eafs))
  stopifnot(length(prs) == nrow(local_anc_counts[[1]]))
  stopifnot(ncol(local_anc_counts[[1]]) == nrow(eafs))
  
 

  EPRS_sum <- local_anc_counts$a1_counts %*% (eafs[,"a1"] * betas_vec) + 
          local_anc_counts$a2_counts %*% (eafs[,"a2"] * betas_vec) + 
          local_anc_counts$a3_counts %*% (eafs[,"a3"] * betas_vec)
  EPRS <- EPRS_sum/length(betas_vec)
  
  varPRS_sum <- local_anc_counts$a1_counts %*% (eafs[,"a1"] * (1 - eafs[,"a1"]) * betas_vec^2) +
                local_anc_counts$a1_counts %*% (eafs[,"a2"] * (1 - eafs[,"a2"]) * betas_vec^2) +
                local_anc_counts$a1_counts %*% (eafs[,"a3"] * (1 - eafs[,"a3"]) * betas_vec^2)
  varPRS <- varPRS_sum/(length(betas_vec)^2)
    
  
  
  rPRS <- prs - EPRS
  
  return(list(lEPRS = EPRS,
              lrPRS = rPRS,
              lvarPRS = varPRS))
  
}
