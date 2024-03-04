
# for n people and one genotype (over one chromosomal copy), sample
# genotype according to allele frequecies of the variant in the ancestral
# population of the variant. 
sample_alleles_one_geno <- function(ancestries_vec, eafs){
  stopifnot(all(ancestries_vec %in% names(eafs)))
  individual_eafs <- eafs[ancestries_vec]
  allele_counts <- rbinom(length(individual_eafs), 1, as.numeric(individual_eafs))
  return(allele_counts)
}
