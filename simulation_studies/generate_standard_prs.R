
# generate standard PRS: weighted sum of alleles 
# (further divided by the number of SNPs used), regardless of
# genetic ancestry.

generate_standard_prs <- function(allele_count_mat, betas_vec){
  stopifnot(length(betas_vec) == ncol(allele_count_mat))
  
  return(drop((allele_count_mat %*% betas_vec)/length(betas_vec)))
  
}