###########################################
## function for simulations
## generate observed PRS and unknown confounders
## this function may return: 
## observed PRS: homogeneous weighting PRS and heterogeneous weighting PRS
## five confounders
## 1. conf-PRS1
## 2. conf-PRS2
## 3. single variant
## 4. two variants
## 5. conf-pc*
## local ancestry counts (for computing local ePRS purposes)
###########################################
###########################################
## load in all the R functions for generating data
###########################################
source("simulate_global_ancestry_pattern.R")
source("sample_local_ancestries.R")
source("sample_local_alleles.R")
source("allele_counts_across_copies.R")
source("generate_standard_prs.R")
source("generate_prs_ancestry_specific_effects.R")
###########################################


simulate_component_data <- function(n, 
                                    ancestral_eaf_for_prs_file, 
                                    ancestral_eaf_for_conf_pc_file,
                                    ancestral_eaf_for_ancestral_confounder_file,
                                    betas_vec_for_standard_prs,
                                    betas_mat_for_prs_ancestry_effects,
                                    admixture_pattern = "realistic_admixture"
){
  
  ## generate global ancestry proportion
  admixed_prop <- simulate_global_ancestry_pattern(n, admixture_pattern)


  ## generate observed PRS -- both homogeneous weighting PRS and heterogeneous weighting PRS
  ## read in the summary statistics and ancestry-specific allele freq data
  prs_eafs <- read.table(ancestral_eaf_for_prs_file, header = TRUE)
  prs_eafs <- prs_eafs[,c("Europe_est", "Africa_est", "Amer_est")]
  colnames(prs_eafs) <- c("a1", "a2", "a3")
  
  ## generate local ancestry 
  local_anc.PRS <- sample_local_ancestries(admixed_prop, p)
  alleles <- sample_local_alleles(local_anc.PRS, prs_eafs)
  ## generate SNP data (count data: 0,1,2)
  allele_counts <- allele_counts_across_copies(alleles)

  ## calculate genetic PCs based on the SNP data that compute observed PRS
  D.PC <- princomp(allele_counts, cor = TRUE)
  
  zdat <- scale(as.matrix(allele_counts))
  
  ## pc score for each samples
  pca.score <- zdat %*% (D.PC$loadings)
  
  ## in simulation 
  ## we consider (1) top 10 PCs; (2) top 20 PCs for method comparison
  Adj.PC10 <- pca.score[, 1:10]
  Adj.PC20 <- pca.score[, 1:20]

  ## generate homogeneous weighting PRS
  standard_prs <- generate_standard_prs(allele_counts, betas_vec_for_standard_prs)
  
  ## generate heterogeneous weighting PRS
  ancestry_effect_prs <- generate_prs_ancestry_specific_effects(local_anc.PRS, 
                                                                alleles,
                                                                betas_mat_for_prs_ancestry_effects)
  
  ## generate unobserved confounding
  ## conf-PRS1
  ## note1: use observed PRS's local ancestry to generate 
  ## note2: use prs_eafs (AF from UKB + ICBP data) to generate genetic data
  ## unknown confounding and observed PRS share same local ancestry information
  ## weighting is generated from N(0, 1)
  
  alleles <- sample_local_alleles(local_anc.PRS, prs_eafs)
  allele_counts <- allele_counts_across_copies(alleles)
  PRS_weighting <- rnorm(p, 0, 1)
  conf_PRS1 <- generate_standard_prs(allele_counts, PRS_weighting)
  
  
  ## generate unobserved confounding
  ## conf-PRS2
  ## based on variants that are common in ancestry "a2" (African)
  ## use MVP summary statistics (AA enriched dataset)
  conf_prs_eafs <- read.table(ancestral_eaf_for_ancestral_confounder_file, header = TRUE)
  conf_prs_effects <- conf_prs_eafs[, "BETA"]
  conf_prs_eafs <- conf_prs_eafs[,c("Europe_est", "Africa_est", "Amer_est")]
  colnames(conf_prs_eafs) <- c("a1", "a2", "a3")
  local_anc.unobsPRS <- sample_local_ancestries(admixed_prop, p)
  alleles <- sample_local_alleles(local_anc.unobsPRS, conf_prs_eafs)
  allele_counts <- allele_counts_across_copies(alleles)
  conf_PRS2 <- generate_standard_prs(allele_counts, conf_prs_effects)
  
  
  ## single variant confounder
  ## generate a confounder variable which is a single variant
  ## randomly select one SNP 
  ## which AF of AA population ranges from 0.45 ~ 0.55
  ind_var <- sample(which(abs(conf_prs_eafs$a2-0.5) < 0.05), 1)
  single_conf_var <- allele_counts[,ind_var]
  
  ## two variants confounder
  ## confounder variables which contains two variants
  ## randomly select two SNPs 
  ## which AF of AA population ranges from 0.45 ~ 0.55
  ind_var <- sample(which(abs(conf_prs_eafs$a2-0.5) < 0.05), 2)
  two_conf_vars <- allele_counts[,ind_var]
  
  
  ## generate conf-pc*
  pc_eafs <- read.table(ancestral_eaf_for_conf_pc_file, header = TRUE)
  pc_eafs <- pc_eafs[,c("Europe_est", "Africa_est", "Amer_est")]
  colnames(pc_eafs) <- c("a1", "a2", "a3")
  p <- nrow(pc_eafs)
  local_anc.pc <- sample_local_ancestries(admixed_prop, p)
  alleles <- sample_local_alleles(local_anc.pc, pc_eafs)
  allele_counts <- allele_counts_across_copies(alleles)
  
  ## simulate the weight from standard normal distribution
  pc_loading <- rnorm(p, 0, 1)
  conf_pc <- generate_standard_prs(allele_counts, pc_loading)
  

  ## output local ancestry
  ## these are all the building blocks for simulation of phenotype-PRS association
  ## also report building blocks needed to compute local ePRS
  ## for local ePRS: local ancestry counts
  ## for each person and each variant, need to report local ancestries
  ## these could be in the form 11, 12, 13, 22, 23, 33 

  local_anc_unphased <- merge_ancestry_alleles_to_unphased(local_anc.PRS)
  local_anc_counts <- merge_ancestry_alleles_to_counts(local_anc.PRS)
  
  return(list(global_ancestries = admixed_prop,
              Data.PC10 = Adj.PC10,
              Data.PC20 = Adj.PC20,
              prs_eafs = prs_eafs,
              standard_prs = standard_prs,
              ancestry_effect_prs = ancestry_effect_prs,
              conf_PRS1 = conf_PRS1,
              conf_PRS2 = conf_PRS2, 
              single_conf_var = single_conf_var, 
              two_conf_vars = two_conf_vars,
              conf_pc = conf_pc,
              local_anc_unphased = local_anc_unphased,
              local_anc_counts = local_anc_counts))
  
}


merge_ancestry_alleles_to_unphased <- function(local_ancestry_alleles){
  vec_local_options <- c("a1a1" = 11,
                         "a1a2" = 12,
                         "a2a1" = 12,
                         "a1a3" = 13,
                         "a3a1" = 12,
                         "a2a2" = 22,
                         "a2a3" = 23,
                         "a3a2" = 23,
                         "a3a3" = 33)
  
  n <- nrow(local_ancestry_alleles[[1]])
  p <- ncol(local_ancestry_alleles[[1]])
  ancestry_alleles_unphased <- matrix(data = NA, nrow = n, ncol = p)
  rownames(ancestry_alleles_unphased) <- rownames(local_ancestry_alleles[[1]])
  
  for (i in 1:p){
    pasted_ancestry_strings <- paste0(local_ancestry_alleles[[1]][,i],
                                      local_ancestry_alleles[[2]][,i])
    ancestry_alleles_unphased[,i] <- vec_local_options[pasted_ancestry_strings]
    
  }
  return(ancestry_alleles_unphased)
}





merge_ancestry_alleles_to_counts <- function(local_ancestry_alleles){
  a1_options <- c("a1a1" = 2,
                  "a1a2" = 1,
                  "a2a1" = 1,
                  "a1a3" = 1,
                  "a3a1" = 1,
                  "a2a2" = 0,
                  "a2a3" = 0,
                  "a3a2" = 0,
                  "a3a3" = 0)
  
  a2_options <- c("a1a1" = 0,
                  "a1a2" = 1,
                  "a2a1" = 1,
                  "a1a3" = 0,
                  "a3a1" = 0,
                  "a2a2" = 2,
                  "a2a3" = 1,
                  "a3a2" = 1,
                  "a3a3" = 0)
  
  a3_options <- c("a1a1" = 0,
                  "a1a2" = 0,
                  "a2a1" = 0,
                  "a1a3" = 1,
                  "a3a1" = 1,
                  "a2a2" = 0,
                  "a2a3" = 1,
                  "a3a2" = 1,
                  "a3a3" = 2)
  
  n <- nrow(local_ancestry_alleles[[1]])
  p <- ncol(local_ancestry_alleles[[1]])
  a1_counts <- a2_counts <- a3_counts <- matrix(data = NA, nrow = n, ncol = p)
  rownames(a1_counts) <- rownames(a2_counts) <- rownames(a3_counts) <- rownames(local_ancestry_alleles[[1]])
  
  for (i in 1:p){
    pasted_ancestry_strings <- paste0(local_ancestry_alleles[[1]][,i],
                                      local_ancestry_alleles[[2]][,i])
    a1_counts[,i] <- a1_options[pasted_ancestry_strings]
    a2_counts[,i] <- a2_options[pasted_ancestry_strings]
    a3_counts[,i] <- a3_options[pasted_ancestry_strings]
    
  }
  return(list(a1_counts = a1_counts,
              a2_counts = a2_counts,
              a3_counts = a3_counts))
}
