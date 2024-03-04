
allele_counts_across_copies <- function(local_alleles){
  stopifnot(all(names(local_alleles) == c("copy_1", "copy_2")))

  local_alleles_copy_1 <- local_alleles$copy_1
  local_alleles_copy_2 <- local_alleles$copy_2
  
  stopifnot(all(dim(local_alleles_copy_1) == dim(local_alleles_copy_2)))
  stopifnot(all(is.element(local_alleles_copy_1, c(0,1))))
  stopifnot(all(is.element(local_alleles_copy_2, c(0,1))))
  
  allele_counts <- local_alleles_copy_1 + local_alleles_copy_2
  
  return(allele_counts)
  
}