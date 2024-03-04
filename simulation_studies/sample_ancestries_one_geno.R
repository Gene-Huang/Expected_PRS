
# A function that samples local ancestries of a specific  variant 
# across a set of people, according to their admixture proportions
# assumes uniform distribution of ancestries (the variant is inherited from 
# ancestry a with probability being the proportion of genome from ancestry a)
# currently assumes 3-way admixture
sample_ancestries_one_geno <- function(global_ancestry_pattern){
  n <- nrow(global_ancestry_pattern)
  ancestry_rand_nums <- runif(n, 0,1)
  sampled_ancestries <- rep("a1", n)
  inds_a2 <- which(ancestry_rand_nums >= global_ancestry_pattern$a1 & 
                     ancestry_rand_nums < global_ancestry_pattern$a1 + global_ancestry_pattern$a2 )
  sampled_ancestries[inds_a2] <- "a2"
  inds_a3 <- which( ancestry_rand_nums >= global_ancestry_pattern$a1 + global_ancestry_pattern$a2 )
  sampled_ancestries[inds_a3] <- "a3"
  
  return(sampled_ancestries)
}
