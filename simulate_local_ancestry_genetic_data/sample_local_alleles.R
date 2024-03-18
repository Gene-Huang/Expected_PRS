########################################################
## R function to simulate SNP data according to the local ancestry pattern
## and corresponding ancestry-specific allele frequency
## tow inputs of this function: (1) local_ancestry_pattern; (2) eafs_by_ancestry
## for the local ancestry pattern input (local_ancestry_pattern)
## contains the information of ancestry patterns across chromosomal copy 1 and copy 2
## for each copy: local ancestries for each of n individuals in each of p variants
## for the ancestry-specific allele frequency input (eafs_by_ancestry)
## provides the information of ancestry-specific allele frequency for each variant 
########################################################
########################################################
## load the R function to generate allele count for each variant according to the local ancestry pattern
########################################################
source("sample_alleles_one_geno.R")
########################################################

sample_local_alleles <- function(local_ancestry_pattern, eafs_by_ancestry){
  
  loc_copy_1 <- local_ancestry_pattern$copy_1
  
  loc_copy_2 <- local_ancestry_pattern$copy_2
  
  stopifnot(ncol(loc_copy_1) == nrow(eafs_by_ancestry))
  stopifnot(all(is.element(as.character(unique(matrix(loc_copy_1))), colnames(eafs_by_ancestry))))
  stopifnot(dim(loc_copy_1) == dim(loc_copy_2))
  
  counts_copy_1 <- counts_copy_2 <- matrix(NA, nrow = nrow(loc_copy_1), ncol = ncol(loc_copy_1))
  
  for (i in 1:ncol(counts_copy_1)){
    
    counts_copy_1[,i] <- sample_alleles_one_geno(loc_copy_1[,i], eafs = eafs_by_ancestry[i,])
    
    counts_copy_2[,i] <- sample_alleles_one_geno(loc_copy_2[,i], eafs = eafs_by_ancestry[i,])
  }
  
  rownames(counts_copy_1) <- rownames(counts_copy_2) <- rownames(loc_copy_1)
  
  return(list(copy_1 = counts_copy_1, copy_2 = counts_copy_2))  
}
