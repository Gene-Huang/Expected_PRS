#########################################################
## Compute_global_ePRS() function
## R function for computing global ePRS
## and the variance of PRS (by global ancestry proportion)
#########################################################
#########################################################
## inputs of this function are:
## input 1: global_ancestry_prop
## global ancestry proportion matrix (n by k + 1)
## the first column is the sample ID information
## n represents number of individuals
## k represents number of ancestries
## the column names of this file should be the same as the example dataset
## which includes: "sample.id", "Africa_est", "SouthernAsian_est", "EasternAsia_est", "Europe_est", "Amer_est", "MiddleEast_est"
## input 2: sumstat
## summary statistics file
## the column names of this file should be the same as the example dataset
## which includes: "rsID", "chr", "hg38_position", "effect_allele", "other_allele", "weight" 
## input 3: ancestry_specific_AF
## ancestry-specific allele frequency file
## the column names of this file should be the same as the example dataset
## which includes: "rsID", "chr", "hg38_position", "REF", "ALT", "Africa_est", "SouthernAsian_est", "EasternAsia_est", "Europe_est", "Amer_est", "MiddleEast_est"
#########################################################
#########################################################
## outputs of "Compute_global_ePRS()" function is a "list"
## the list contains two outputs
## first output: global ePRS value for each individual
## second output: variance of PRS for each individual
#########################################################

#########################################################
## load R packages
#########################################################
library(dplyr)

library(data.table)

#########################################################
## load R function
#########################################################
## R function for computing global ePRS
#########################################################
global_ePRS <- function(gaProp, sumstat, ancestry_specific_af){
  
  ## step 1
  ## sumstat*ancestry-specific allele frequency
  ## sum across all variants
  ## times 2
  weighted_af <- 2*colSums(ancestry_specific_af*sumstat)
  
  ## step2: sum over ancestry weighted by global ancestry proportion
  final_global_ePRS <- gaProp %*% weighted_af 
  
  return(final_global_ePRS)
  
}

#########################################################
## R function for computing variance of PRS
## by global ancestry proportion
#########################################################
variance_PRS_gaProp <- function(gaProp, sumstat, ancestry_specific_af){
  
  ## step 1
  ## (sumstat)^2*ancestry-specific allele frequency
  ## sum across all variants
  ## times 2
  weighted_af <- colSums((2*ancestry_specific_af*(1-ancestry_specific_af))*(sumstat^2))
  
  ## step2: sum over ancestry weighted by global ancestry proportion
  final_variance_PRS_gaProp <- gaProp %*% weighted_af 
  
  return(final_variance_PRS_gaProp)
  
}


#########################################################
## start the "Compute_global_ePRS()" function
## R function for computing global ePRS
## and the variance of PRS (by global ancestry proportion)
#########################################################

Compute_global_ePRS <- function(global_ancestry_prop, sumstat, ancestry_specific_AF){

  #########################################
  ## for sumstat file
  #########################################
  ## step1: 
  ## delete duplicated variants
  ## by rsID
  #########################################
  sumstat <- sumstat[!(duplicated(sumstat$rsID) | duplicated(sumstat$rsID, fromLast = TRUE)), ]
  
  #########################################
  ## step2: 
  ## create two columns in sumstat file
  ## in order to match to ancestry-specific AF file
  ## 1st variant ID: chromosome:hg38_position:effect_allele:other_allele
  ## 2nd variant ID: chromosome:hg38_position:other_allele:effect_allele
  #########################################
  sumstat$chr_pos_effect_other <- paste0(sumstat$chr, ":",
                                         sumstat$hg38_position, ":",
                                         sumstat$effect_allele, ":",
                                         sumstat$other_allele)
  
  
  sumstat$chr_pos_other_effect <- paste0(sumstat$chr, ":",
                                         sumstat$hg38_position, ":",
                                         sumstat$other_allele, ":",
                                         sumstat$effect_allele)
  
  
  #########################################
  ## step3:
  ## delete duplicated variants
  ## just for double-checking purposes
  ## by chr:pos:effect:other
  #########################################
  sumstat <- sumstat[!(duplicated(sumstat$chr_pos_effect_other) | duplicated(sumstat$chr_pos_effect_other, fromLast = TRUE)), ]
  
  
  #########################################
  ## for ancestry-specific AF file
  #########################################
  ## step1: 
  ## delete duplicated variants
  ## by rsID
  #########################################
  ancestry_specific_AF <- ancestry_specific_AF[!(duplicated(ancestry_specific_AF$rsID) | duplicated(ancestry_specific_AF$rsID, fromLast = TRUE)), ]
  
  #########################################
  ## step2
  ## create one column in ancestry-specific AF file
  ## in order to match to sumstat file
  ## the variant ID: chromosome:hg38_position:alt_allele:ref_allele
  #########################################
  ancestry_specific_AF$chr_pos_alt_ref <- paste0(ancestry_specific_AF$chr, ":",
                                                 ancestry_specific_AF$hg38_position, ":",
                                                 ancestry_specific_AF$ALT, ":",
                                                 ancestry_specific_AF$REF)

  #########################################
  ## step3:
  ## delete duplicated variants
  ## just for double-checking purposes
  ## by chr:pos:alt:ref
  #########################################
  ancestry_specific_AF <- ancestry_specific_AF[!(duplicated(ancestry_specific_AF$chr_pos_alt_ref) | duplicated(ancestry_specific_AF$chr_pos_alt_ref, fromLast = TRUE)), ]
  
  
  #########################################
  ## Note: important step!!!
  ## before computing ePRS
  ## we need to check that
  ## the effect allele used in GWAS (the weight was computed based on)
  ## matches the allele that the allele frequency was calculated based on in TOPMed
  #########################################
  #########################################
  ## start to do matching
  ## check that the effect allele in sumstat file matched to ALT or matched to REF in ancestry-specific AF file
  ## Note: in the ancestry-specific AF file
  ## the allele freq was computed based the "ALT" allele
  #########################################
  ## effect allele is ALT
  AF_match_effect <- ancestry_specific_AF[which(ancestry_specific_AF$chr_pos_alt_ref %in% sumstat$chr_pos_effect_other),]
  AF_match_effect$match_info <- "alt"
  
  ## effect allele is REF
  AF_match_other <- ancestry_specific_AF[which(ancestry_specific_AF$chr_pos_alt_ref %in% sumstat$chr_pos_other_effect),]
  AF_match_other$match_info <- "ref"
  
  ## combine these two datasets together
  final_matched_AF <- rbind(AF_match_effect, AF_match_other)
  
  
  ## merge summary stat file and allele freq file
  merge_data <- left_join(final_matched_AF, sumstat, by = "rsID")

  ## select the variables we may use to compute ePRS
  merge_data <- merge_data[,c("rsID", "chr_pos_effect_other", "chr.x", "hg38_position.x",
                              "REF", "ALT", "effect_allele", "other_allele", "match_info", 
                              "weight",
                              "Africa_est", "SouthernAsian_est", "EasternAsia_est", "Europe_est", "Amer_est", "MiddleEast_est")]
  
  
  ## rename the variable names (not necessary)
  merge_data <- rename(merge_data, 
                       "chr" = "chr.x",
                       "hg38_position" = "hg38_position.x")
  
  
  #########################################
  ## start to compute global ePRS
  #########################################
  #########################################
  ## step1: 
  ## prepare the input file of ancestry-specific allele frequency
  #########################################
  ## select 6 ancestries' AF
  input_AF <- merge_data[, c("Africa_est", "SouthernAsian_est", "EasternAsia_est", "Europe_est", "Amer_est", "MiddleEast_est")]
  
  ## because the AF computed from TOPMed was based on "ALT" allele
  ## so if the effect allele from the summary statistics file was matched to "Ref"
  ## the ancestry allele freq should be recalculated as: 1-AF
  input_AF[which(merge_data$match_info == "ref"),] <- 1 - input_AF[which(merge_data$match_info == "ref"),]
  
  #########################################
  ## step2: 
  ## prepare the input file of weight (summary statistics)
  #########################################
  input_sumstat <- merge_data$weight
  
  #########################################
  ## step3: 
  ## prepare the input file of global ancestry proportion
  #########################################
  sample_ID <- global_ancestry_prop$sample.id
  
  global_ancestry_prop <- global_ancestry_prop[, -"sample.id"]
  
  ## make sure the ancestry proportion can match to the ancestry-specific AF file 
  ## the order of the colnames of these two datasets should be the same
  global_ancestry_prop <- global_ancestry_prop[,c("Africa_est", "SouthernAsian_est", "EasternAsia_est", "Europe_est", "Amer_est", "MiddleEast_est")]
  
  input_gaProp <- matrix(unlist(global_ancestry_prop), nrow(global_ancestry_prop), ncol(global_ancestry_prop))
  
  #########################################
  ## compute global ePRS
  ## use global_ePRS() function
  #########################################
  gePRS <- global_ePRS(gaProp = input_gaProp,
                       sumstat = input_sumstat,
                       ancestry_specific_af = input_AF)
  
  
  gePRS <- data.frame(sample.id = sample_ID, global_ePRS = gePRS)
  
  #########################################
  ## compute variance of PRS
  ## use variance_PRS_gaProp() function
  #########################################
  var_global_ePRS <- variance_PRS_gaProp(gaProp = input_gaProp,
                                         sumstat = input_sumstat,
                                         ancestry_specific_af = input_AF)
  
  var_global_ePRS <- data.frame(sample.id = sample_ID, variance_PRS = var_global_ePRS)
  
  
  ## output
  return(list(global_ePRS = gePRS, variance_PRS = var_global_ePRS))
  
}



