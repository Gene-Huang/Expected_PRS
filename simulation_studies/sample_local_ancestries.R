
# sample local ancestry for all people represented by a 
# matrix global_ancestry_pattern, in p variants.

source("sample_ancestries_one_geno.R")

sample_local_ancestries <- function(global_ancestry_pattern, p){
  
  sampled_ancestries_copy_1 <- replicate(p, sample_ancestries_one_geno((global_ancestry_pattern)))
  sampled_ancestries_copy_2 <- replicate(p, sample_ancestries_one_geno((global_ancestry_pattern)))
  rownames(sampled_ancestries_copy_1) <- rownames(sampled_ancestries_copy_2) <- rownames(global_ancestry_pattern)
  
    return(list(copy_1 = sampled_ancestries_copy_1, 
                copy_2 = sampled_ancestries_copy_2))
}


