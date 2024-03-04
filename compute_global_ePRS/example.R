##################################################
## examples for computing global ePRS
## and the variance of PRS
##################################################
## load R package
library(data.table)

##################################################
## prepare the input files for using Compute_global_ePRS() function
## Note: all the column names of the input files
## should be the same as the example files
##################################################
##################################################
## input-1
## global ancestry proportion
## this dataset is simulated data for testing purposes
## includes 100 individuals
## global ancestry proportion of 6 ancestries
##################################################
input1 <- fread("gaProp_example.csv")

head(input1)

##################################################
## input-2
## summary statistics file
## subset of BMI GWAS summary statistics
##################################################
input2 <- fread("sumstat_example.csv")

head(input2)

##################################################
## input-3
## ancestry-specific allele frequency file
## subset of the ancestry-specific allele frequency computed using TOPMed individuals
## note: the allele frequency was computed based the "ALT" allele in TOPMed
##################################################
input3 <- fread("ancestry_specific_AF_by_GAFA_example.csv")

head(input3)

##################################################
## load in the R function for computing global ePRS and variance of PRS
##################################################
source("compute_global_ePRS.R")

Output <- Compute_global_ePRS(global_ancestry_prop = input1, 
                              sumstat = input2, 
                              ancestry_specific_AF = input3)


## global ePRS of each individual
Output$global_ePRS


## variance of PRS of each individual
Output$variance_PRS


