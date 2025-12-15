##########################################
## 2025.12.15
## step1: prepare variant list
##########################################

## load R packages
library(SeqArray)
library(SeqVarTools)
library(gdsfmt)
library(data.table)

## input data
## from 1000 Genome Phase 3 reference panel
## https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

gds_chr2 <- "1000G_chr2.gds"
panel_fp <- "integrated_call_samples_v3.20130502.ALL.panel.txt"

target_n <- 100000
min_call_rate <- 0.99
maf_min <- 0.01
drop_palindromic <- TRUE

rng_seed <- 12152025       
enforce_phase3_n <- FALSE
phase3_n_expected <- 2504

## load GDS file
if (!file.exists(gds_chr2)) stop("GDS not found: ", gds_chr2)
if (!file.exists(panel_fp)) stop("Panel not found: ", panel_fp)
if (!is.null(rng_seed)) set.seed(as.integer(rng_seed))

g <- seqOpen(gds_chr2)
on.exit(try(seqClose(g), silent = TRUE), add = TRUE)

samples <- seqGetData(g, "sample.id")
if (enforce_phase3_n) stopifnot(length(samples) == phase3_n_expected)

panel <- fread(panel_fp)
setnames(panel, c("sample","pop","super_pop","gender"), c("sample","pop","super","sex"))
panel <- panel[match(samples, sample)]
if (anyNA(panel$sample)) stop("Some GDS samples were not found in the panel file.")

grp <- list(
  AFR = which(panel$super == "AFR"),
  EUR = which(panel$super == "EUR"),
  EAS = which(panel$super == "EAS")
)
stopifnot(all(vapply(grp, length, 1L) > 0))

## focus on chr2
## rsID 
## biallelic SNPs (A/C/G/T)
seqResetFilter(g, verbose = FALSE)
chrom_all <- seqGetData(g, "chromosome")
idx_chr2 <- which(chrom_all %in% c("2","chr2",2L))
if (!length(idx_chr2)) stop("No chr2 variants found.")

vid_all <- seqGetData(g, "variant.id")
seqSetFilter(g, variant.sel = vid_all[idx_chr2], verbose = FALSE)

## "REF,ALT"
allele <- seqGetData(g, "allele")   
alleles <- strsplit(allele, ",", fixed = TRUE)
rsid <- seqGetData(g, "annotation/id")
has_rsid <- !is.na(rsid) & nzchar(rsid)
biallel <- vapply(alleles, function(a) length(a) == 2L, TRUE)
is_base <- vapply(alleles, function(a) all(nchar(a) == 1L & a %in% c("A","C","G","T")), TRUE)
keep0 <- biallel & is_base & has_rsid

if (!any(keep0)) stop("No biallelic rsID SNPs on chr2 after filtering.")
cur_vid <- seqGetData(g, "variant.id")
seqSetFilter(g, variant.sel = cur_vid[keep0], verbose = FALSE)

## drop palindromic SNPs (A/T and C/G)
if (drop_palindromic) {
  refalt <- seqGetData(g, "allele")
  ra <- tstrsplit(refalt, ",", fixed = TRUE)
  ref <- ra[[1]]; alt <- ra[[2]]
  is_pal <- (ref == "A" & alt == "T") | (ref == "T" & alt == "A") |
    (ref == "C" & alt == "G") | (ref == "G" & alt == "C")
  n_drop <- sum(is_pal, na.rm = TRUE)
  if (n_drop > 0L) {
    cur_vid <- seqGetData(g, "variant.id")
    seqSetFilter(g, variant.sel = cur_vid[!is_pal], verbose = FALSE)
    message(sprintf("Dropped %d palindromic SNPs.", n_drop))
  } else {
    message("No palindromic SNPs found.")
  }
}

## get exact variant.id 
cur_vid_at_stats <- seqGetData(g, "variant.id")

## call rate and AF of ATL allele
.per_group_stats <- function(idx) {
  seqSetFilter(g, sample.id = samples[idx], verbose = FALSE)
  list(
    callrate = 1 - missingGenotypeRate(g),
    af = alleleFrequency(g, n = 1)  
  )
}

st_AFR <- .per_group_stats(grp$AFR)
st_EUR <- .per_group_stats(grp$EUR)
st_EAS <- .per_group_stats(grp$EAS)


st_ALL <- (function() {
  seqSetFilter(g, sample.id = samples, verbose = FALSE)
  list(
    callrate = 1 - missingGenotypeRate(g),
    af = alleleFrequency(g, n = 1)  
  )
})()


seqSetFilter(g, sample.id = NULL, verbose = FALSE)

## variants should
## (1) present in three populations
## (2) pass call rate screening
## (3) pass MAF screening 

present_all3 <-
  (st_AFR$callrate >= min_call_rate) &
  (st_EUR$callrate >= min_call_rate) &
  (st_EAS$callrate >= min_call_rate)

maf_AFR <- pmin(st_AFR$af, 1 - st_AFR$af)
maf_EUR <- pmin(st_EUR$af, 1 - st_EUR$af)
maf_EAS <- pmin(st_EAS$af, 1 - st_EAS$af)

pass_mask <- present_all3 &
  (maf_AFR >= maf_min) &
  (maf_EUR >= maf_min) &
  (maf_EAS >= maf_min)

if (!any(pass_mask)) stop("No variants pass AFR/EUR/EAS call-rate & MAF filters.")

## randomly select variants for further analysis
pass_vid <- cur_vid_at_stats[pass_mask]
if (length(pass_vid) < target_n)
  stop(sprintf("Only %d variants meet criteria; need %d.", length(pass_vid), target_n))

pick_vid <- sort(sample(pass_vid, target_n, replace = FALSE))  
seqSetFilter(g, variant.sel = pick_vid, verbose = FALSE)

## prepare final data (variant information)
pos <- seqGetData(g, "position")
refalt <- seqGetData(g, "allele")
rsid_pick <- seqGetData(g, "annotation/id")
chrom_col <- rep("2", length(pos))

ra <- tstrsplit(refalt, ",", fixed = TRUE)
ref <- ra[[1]]
alt <- ra[[2]]

idx_in_stats <- match(pick_vid, cur_vid_at_stats)

Eur_af_alt <- st_EUR$af[idx_in_stats]
Afr_af_alt <- st_AFR$af[idx_in_stats]
EAS_af_alt <- st_EAS$af[idx_in_stats]
ALL_af_alt <- st_ALL$af[idx_in_stats]

ALL_maf <- pmin(ALL_af_alt, 1 - ALL_af_alt)
ALL_minor <- ifelse(ALL_af_alt < 0.5, alt,
                    ifelse(ALL_af_alt > 0.5, ref, "TIE"))

out_final <- data.table(
  chrom = chrom_col,
  pos = pos,
  rsid = rsid_pick,
  ref = ref,
  alt = alt,
  ALL_minor = ALL_minor,
  ALL_maf = ALL_maf,
  Eur_af_alt = Eur_af_alt,
  Afr_af_alt = Afr_af_alt,
  EAS_af_alt = EAS_af_alt
)

head(out_final)



