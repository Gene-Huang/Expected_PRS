########################################################
## R function to generate allele count for each variants
## according to the local ancestry pattern and the ancestry-specific allele frequency
## for n individual and one genotype (over one chromosomal copy)
## sample allele count according to allele frequencies of the variant in the ancestral population of the variant
########################################################

sample_alleles_one_geno <- function(ancestries_vec, eafs){
  
  stopifnot(all(ancestries_vec %in% names(eafs)))
  
  individual_eafs <- eafs[ancestries_vec]
  
  ## for each copy, for each variant
  ## simulate the allele count from Bernoulli distribution
  ## with the probability equal to ancestry-specific allele frequency
  allele_counts <- rbinom(length(individual_eafs), 1, as.numeric(individual_eafs))
  
  return(allele_counts)
}
