
# generate EPRS for the standard PRS: weighted sum of alleles 
# (further divided by the number of SNPs used), regardless of
# genetic ancestry.

generate_global_EPRS_rPRS <- function(prs, 
                                      eafs, 
                                      betas_vec,
                                      admixed_prop){
  stopifnot(length(betas_vec) == nrow(eafs))
  stopifnot(length(prs) == nrow(admixed_prop))
  stopifnot(ncol(admixed_prop) == ncol(eafs))
  
  Eprs_full_ancestry <- varprs_full_ancestry <- rep(NA, ncol(admixed_prop))
  names(Eprs_full_ancestry) <- names(varprs_full_ancestry) <- colnames(admixed_prop)
  
  # compute EPRS for each ancestry.
  for (i in 1:length(Eprs_full_ancestry)){
    cur_ancestry <- names(Eprs_full_ancestry)[i]
    Eprs_full_ancestry[i] <- sum(2*eafs[,cur_ancestry]*betas_vec)/length(betas_vec)
    varprs_full_ancestry[i] <- sum(2*eafs[,cur_ancestry]*(1-eafs[,cur_ancestry])*betas_vec^2)/(length(betas_vec)^2)
  }
  
  # compute expected global PRS per person, by (global) proportion ancestry
  EPRS <- as.matrix(admixed_prop) %*% Eprs_full_ancestry
  varPRS <- as.matrix(admixed_prop) %*% varprs_full_ancestry
  rPRS <- prs - EPRS
  
  return(list(gEPRS = EPRS,
              grPRS = rPRS,
              gvarEPRS = varPRS))
  
}
