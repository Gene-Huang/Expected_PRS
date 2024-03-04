
# generate PRS with ancestry specific effects: weighted sum of alleles 
# (further divided by the number of SNPs used), with weights
# potentially dependent on the ancestry of the genetic variant.

generate_prs_ancestry_specific_effects <- function(local_anc_list,
                                                   allele_count_list,
                                                   betas_mat){
  stopifnot(length(local_anc_list) == 2)
  stopifnot(length(allele_count_list) == 2)
  
  p <- nrow(betas_mat)
  n <- nrow(local_anc_list[[1]])
  stopifnot(p == ncol(local_anc_list[[1]]))
  stopifnot(p == ncol(allele_count_list[[1]]))
  
  # sum contribution of alleles in "copy 1" of the variants
  prs_copy_1 <- prs_copy_2 <- rep(0, n)
  for (i in 1:p){
    prs_copy_1 <- prs_copy_1 + 
                  allele_count_list[[1]][,i]*betas_mat[i,local_anc_list[[1]][,i]]
    prs_copy_2 <- prs_copy_2 + 
      allele_count_list[[2]][,i]*betas_mat[i,local_anc_list[[2]][,i]]
  }
  
  prs <- (prs_copy_1 + prs_copy_2)/p
  
  return(prs)
  
}

