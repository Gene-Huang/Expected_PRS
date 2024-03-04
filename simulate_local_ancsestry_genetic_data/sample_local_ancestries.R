###############################################
## sample local ancestry for all individuals
## across all p variants
## based on global ancestry proportion file
###############################################
###############################################
## load the function for generating local ancestry of a sible variant
###############################################
source("sample_ancestries_one_geno.R")
###############################################

sample_local_ancestries <- function(global_ancestry_pattern, p){
  
  ## sample the first copy
  sampled_ancestries_copy_1 <- replicate(p, sample_ancestries_one_geno((global_ancestry_pattern)))
  
  ## sample the second copy
  sampled_ancestries_copy_2 <- replicate(p, sample_ancestries_one_geno((global_ancestry_pattern)))
  
  rownames(sampled_ancestries_copy_1) <- rownames(sampled_ancestries_copy_2) <- rownames(global_ancestry_pattern)
  
  return(list(copy_1 = sampled_ancestries_copy_1,
              copy_2 = sampled_ancestries_copy_2))
  
}


