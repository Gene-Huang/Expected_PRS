##########################################
## 2025.12.15
## step2: PRS-outcome association analysis
## continuous traits
##########################################

## load R packages
library(SeqArray)
library(SeqVarTools)
library(gdsfmt)
library(data.table)


## 1000 genome phase 3
## 2504 individuals
## AFR AMR EAS EUR SAS 
## 661 347 504 503 489 

## load in ePRS function
source("generate_global_EPRS_rPRS.R")

source("generate_local_EPRS_rPRS.R")


## input datasets
## "chr2_100k_random_variants_from_1000G.csv": output from step 1

gds_fp <- "1000G_chr2.gds"  
panel_fp <- "integrated_call_samples_v3.20130502.ALL.panel.txt"
variants_fp <- "chr2_100k_random_variants_from_1000G.csv" 


## set simulation parameters
## number of individuals
N_indiv <- 20000

## number of SNPs to simulate 
m_per_indiv <- 1000

## expected tract length (in SNPs) for geometric dist.
mean_block <- 10


## for checking purposes
strict_sanity_checks <- TRUE

## Dirichlet distribution
## for generating global ancestry proportion
rdirichlet1 <- function(alpha) { x <- rgamma(length(alpha), alpha, 1); x/sum(x) }

draw_props <- function(
    scenario_weights = c(`EUR-major`=0.50, `AFR-major`=0.20, `EAS-major`=0.15, `Tri-hybrid`=0.15),
    floors = c(EUR=0.01, AFR=0.01, EAS=0.01),
    concentration = 12
) {
  scen <- sample(names(scenario_weights), 1, prob = scenario_weights)
  mu <- switch(scen,
               "EUR-major" = c(EUR=0.70, AFR=0.15, EAS=0.15),
               "AFR-major" = c(EUR=0.20, AFR=0.70, EAS=0.10),
               "EAS-major" = c(EUR=0.05, AFR=0.10, EAS=0.85),
               "Tri-hybrid" = c(EUR=0.40, AFR=0.30, EAS=0.30))
  p <- rdirichlet1(concentration * mu)
  p <- pmax(p, floors); p <- p / sum(p)
  setNames(p, c("EUR", "AFR", "EAS"))
}

## extract haplotype 
.extract_haps <- function(g, samples, sample_idx, m_expected, strict = TRUE) {
  
  seqSetFilter(g, sample.id = samples[sample_idx], verbose = FALSE)
  GT <- seqGetData(g, "genotype")
  dm <- dim(GT)
  if (length(dm) != 3) {
    stop(sprintf("Genotype array is not 3D. dim=%s", paste(dm, collapse="x")))
  }
  
  n_samp <- length(sample_idx)
  ax_sample <- which(dm == n_samp)
  ax_variant <- which(dm == m_expected)
  
  if (length(ax_sample) != 1L || length(ax_variant) != 1L) {
    stop(sprintf("Cannot infer sample/variant axes from dim=%s (n_samp=%d, m_expected=%d).",
                 paste(dm, collapse="x"), n_samp, m_expected))
  }
  ax_ploidy <- setdiff(1:3, c(ax_sample, ax_variant))
  if (length(ax_ploidy) != 1L) {
    stop(sprintf("Cannot resolve ploidy axis from dim=%s.", paste(dm, collapse="x")))
  }
  if (dm[ax_ploidy] != 2L) {
    stop(sprintf("Ploidy != 2 (found %d) — not supported for autosomes.", dm[ax_ploidy]))
  }
  
  GT_svp <- aperm(GT, c(ax_sample, ax_variant, ax_ploidy))
  
  GT_svp[GT_svp > 1L] <- 1L
  
  n <- dim(GT_svp)[1]; m <- dim(GT_svp)[2]
  if (!missing(m_expected)) stopifnot(m == m_expected)

  H <- matrix(0L, nrow = 2L*n, ncol = m)
  H[seq(1, 2*n, by=2), ] <- GT_svp[,,1]
  H[seq(2, 2*n, by=2), ] <- GT_svp[,,2]
  
  if (strict) {
    ds   <- altDosage(g)
    stopifnot(nrow(ds) == n, ncol(ds) == m)
    Dhat <- H[seq(1, 2*n, by=2), ] + H[seq(2, 2*n, by=2), ]
    if (max(abs(Dhat - ds), na.rm = TRUE) != 0) {
      stop("ALT dosage mismatch between hap sum and altDosage().")
    }
  }
  H
}

## simulate local ancestry state
draw_tracts <- function(m, props, mean_block=50L) {
  labs <- names(props); p <- as.numeric(props)/sum(props)
  labs_out <- character(m); i <- 1L
  while (i <= m) {
    a <- sample(labs, 1L, prob = p)
    L <- rgeom(1L, prob = 1/mean_block) + 1L
    j <- min(m, i + L - 1L)
    labs_out[i:j] <- a
    i <- j + 1L
  }
  labs_out
}

copy_hap <- function(anc_vec, H_EUR, H_AFR, H_EAS) {
  m <- length(anc_vec)
  r <- rle(anc_vec)
  hap <- integer(m)
  idx <- 1L
  for (k in seq_along(r$lengths)) {
    a <- r$values[k]; L <- r$lengths[k]; j <- idx + L - 1L
    donor_row <- switch(a,
                        EUR = sample.int(nrow(H_EUR), 1L),
                        AFR = sample.int(nrow(H_AFR), 1L),
                        EAS = sample.int(nrow(H_EAS), 1L))
    hap[idx:j] <- switch(a,
                         EUR = H_EUR[donor_row, idx:j],
                         AFR = H_AFR[donor_row, idx:j],
                         EAS = H_EAS[donor_row, idx:j])
    idx <- j + 1L
  }
  hap
}

simulate_one <- function(H_EUR, H_AFR, H_EAS, props, mean_block=50L) {
  m <- ncol(H_EUR)
  anc1 <- draw_tracts(m, props, mean_block)
  anc2 <- draw_tracts(m, props, mean_block)
  hap1 <- copy_hap(anc1, H_EUR, H_AFR, H_EAS)
  hap2 <- copy_hap(anc2, H_EUR, H_AFR, H_EAS)
  list(
    alt_count_sum  = hap1 + hap2,  # 0/1/2 ALT copies
    alt_count_hap1 = hap1,         # 0/1
    alt_count_hap2 = hap2,         # 0/1
    la_hap1        = anc1,         # "EUR"/"AFR"/"EAS"
    la_hap2        = anc2,
    props          = props
  )
}


## start simulation
stopifnot(file.exists(gds_fp), file.exists(panel_fp), file.exists(variants_fp))
g <- seqOpen(gds_fp); on.exit(try(seqClose(g), silent=TRUE), add=TRUE)

samples <- seqGetData(g, "sample.id")
panel <- fread(panel_fp)
setnames(panel, c("sample","pop","super_pop","gender"), c("sample","pop","super","sex"))
panel <- panel[match(samples, sample)]
stopifnot(!anyNA(panel$sample))

grp <- list(
  EUR = which(panel$super == "EUR"),
  AFR = which(panel$super == "AFR"),
  EAS = which(panel$super == "EAS")
  )
stopifnot(all(vapply(grp, length, 1L) > 0))
  
## load in variant list 
dt10k <- fread(variants_fp)
req_cols <- c("rsid","chrom","pos","ref","alt","Eur_af_alt","Afr_af_alt","EAS_af_alt")
if (!all(req_cols %in% names(dt10k))) {
  stop(sprintf("variants_fp must contain columns: %s", paste(req_cols, collapse=", ")))
  }
stopifnot(all(dt10k$chrom == "2"))
  
seqResetFilter(g, verbose=FALSE)
chrom_all <- seqGetData(g, "chromosome")
vid_all <- seqGetData(g, "variant.id")
idx_chr2 <- which(chrom_all %in% c("2","chr2",2L))
stopifnot(length(idx_chr2) > 0)
seqSetFilter(g, variant.sel = vid_all[idx_chr2], verbose=FALSE)
  

allele_str <- seqGetData(g, "allele")
aa <- strsplit(allele_str, ",", fixed=TRUE)
is_bi <- vapply(aa, function(x) length(x)==2L, TRUE)
is_ATCG <- vapply(aa, function(x) all(nchar(x)==1L & x %in% c("A","C","G","T")), TRUE)
keep_snp <- which(is_bi & is_ATCG)
vid_chr2 <- seqGetData(g, "variant.id")
seqSetFilter(g, variant.sel = vid_chr2[keep_snp], verbose=FALSE)

g_rsid <- seqGetData(g, "annotation/id")
g_vid <- seqGetData(g, "variant.id")
rs_dt <- data.table(variant.id = g_vid, rsid = g_rsid)
rs_dt <- rs_dt[!is.na(rsid) & nzchar(rsid)]
rs_dt <- rs_dt[!duplicated(rsid)]
setkey(rs_dt, rsid)
  
cand <- rs_dt[J(dt10k$rsid), nomatch=0L] 
if (nrow(cand) < m_per_indiv)
  stop(sprintf("Only %d candidate SNPs available after mapping; need %d.", nrow(cand), m_per_indiv))
sel_vid <- sample(cand$variant.id, m_per_indiv, replace = FALSE)
seqSetFilter(g, variant.sel = sel_vid, verbose=FALSE)
ord <- order(seqGetData(g, "position"))
vid_sorted <- seqGetData(g, "variant.id")[ord]
seqSetFilter(g, variant.sel = vid_sorted, verbose = FALSE)
  

pos_cur <- seqGetData(g, "position")
rsid_cur <- seqGetData(g, "annotation/id")
refalt <- seqGetData(g, "allele")
ra <- tstrsplit(refalt, ",", fixed=TRUE)
sel_coords <- data.table(
  chrom = rep("2", length(pos_cur)),
  pos = pos_cur,
  rsid = rsid_cur,
  ref = ra[[1]],
  alt = ra[[2]]
  )
m_now <- nrow(sel_coords)
stopifnot(m_now == m_per_indiv)
message(sprintf("Using %d SNPs after mapping & SNP enforcement.", m_now))
  

## extract ancestry-specific haplotypes (ALT-coded)
H_EUR <- .extract_haps(g, samples, grp$EUR, m_expected = m_now, strict = strict_sanity_checks)
H_AFR <- .extract_haps(g, samples, grp$AFR, m_expected = m_now, strict = strict_sanity_checks)
H_EAS <- .extract_haps(g, samples, grp$EAS, m_expected = m_now, strict = strict_sanity_checks)
seqSetFilter(g, sample.id = NULL, verbose = FALSE)
stopifnot(ncol(H_EUR)==m_now, ncol(H_AFR)==m_now, ncol(H_EAS)==m_now)
  
  
  
## Close the GDS early 
try(seqClose(g), silent = TRUE)
g <- NULL

props_list <- replicate(N_indiv, draw_props(), simplify = FALSE)

res <- lapply(props_list, function(p)
  simulate_one(H_EUR, H_AFR, H_EAS, props = p, mean_block = mean_block))

  
ALT_counts_sum <- do.call(rbind, lapply(res, `[[`, "alt_count_sum"))  
ALT_counts_hap1 <- do.call(rbind, lapply(res, `[[`, "alt_count_hap1")) 
ALT_counts_hap2 <- do.call(rbind, lapply(res, `[[`, "alt_count_hap2")) 
LA_hap1_chr <- do.call(rbind, lapply(res, `[[`, "la_hap1"))    
LA_hap2_chr <- do.call(rbind, lapply(res, `[[`, "la_hap2"))
gaProp <- do.call(rbind, lapply(res, `[[`, "props"))    
  

af_ref <- merge(sel_coords[, .(rsid)],
                dt10k[, .(rsid, Eur_af_alt, Afr_af_alt, EAS_af_alt)],
                by = "rsid", all.x = TRUE, sort = FALSE)
if (anyNA(af_ref$Eur_af_alt) || anyNA(af_ref$Afr_af_alt) || anyNA(af_ref$EAS_af_alt)) {
  
  stop("Missing ancestry-specific ALT AFs for some selected SNPs in variants_fp.")
}


  
## prepare variant information
val_dt <- data.table(
  rsid = sel_coords$rsid,
  chrom = sel_coords$chrom,
  pos = sel_coords$pos,
  ref = sel_coords$ref,
  alt = sel_coords$alt,
  EUR_af_alt = af_ref$Eur_af_alt,
  AFR_af_alt = af_ref$Afr_af_alt,
  EAS_af_alt = af_ref$EAS_af_alt
  )
  

test_af <- val_dt[, c("EUR_af_alt", "AFR_af_alt", "EAS_af_alt")]

## get variants that are afr-enriched (allele freq)
afr_enrich_index <- apply(test_af, 1, function(input) ifelse(input[2] >= 0.8 &
                                                               input[1] < 0.5 &
                                                               input[3] < 0.5, 1, 0))
  
val_dt$afr_enrich_index <- afr_enrich_index
  
val_dt$af_no_enrich_index <- 1- val_dt$afr_enrich_index


  
######################################
## generate phenotype Y
######################################
######################################
## set-up
######################################
  
## target SNP-heritability 
h2 <- 0.50
#h2 <- 0.20
#h2 <- 0.80
  
## proportion causal SNPs
prop_causal <- 0.01 
  
## number of variant in analysis
m <- nrow(val_dt)
  
## number of causal variants
M <- m*prop_causal
  
## number of sample size
N <- nrow(ALT_counts_sum)
  
exact_h2 <- TRUE

## strength for confounding: 0.5, 1.0, 1.5
alpha <- 0.5
#alpha <- 1.0
#alpha <- 1.5

## select causal variant
causal_idx <- sample(which(val_dt$af_no_enrich_index == 1), M, replace = "FALSE")   

p_eur_alt <- val_dt$EUR_af_alt

## generate effect size for causal variants
beta_norm <- rnorm(M, mean = 0, sd = sqrt(h2 / M))
denom <- sqrt(2 * p_eur_alt[causal_idx] * (1 - p_eur_alt[causal_idx]))
b_causal <- beta_norm / denom
b <- numeric(m)
b[causal_idx] <- b_causal
val_dt$beta <- b
  
## compute oracle PRS 
gval <- as.numeric(ALT_counts_sum %*% b)
if (exact_h2) {
  sf <- sqrt(h2 / var(gval))
  b <- b * sf
  gval <- gval * sf
  }
  
## compute phenotype
eps <- rnorm(N, 0, sqrt(1 - h2))

Y <- gval + eps
  
## select single variant confounding
conf_index <- sample(which(val_dt$afr_enrich_index == 1 & val_dt$af_no_enrich_index == 0), 1, replace = "FALSE")

conf <- ALT_counts_sum[, conf_index]

## add confounding
Y <- Y + alpha*(conf/sd(conf))
  
  
######################################
## run GWAS using training data
######################################

n_discovery <- 10000
n_target <- 10000
  
## PCs used for GWAS (training dataset only)
use_pcs <- TRUE

n_pcs <- 5
N <- N_indiv
  
## prepare training and validation set
perm <- sample.int(N, N, replace = FALSE)
idx_dis <- perm[seq_len(n_discovery)]
idx_tgt <- perm[(n_discovery + 1):(n_discovery + n_target)]
G_dis <- ALT_counts_sum[idx_dis, , drop = FALSE]
G_tgt <- ALT_counts_sum[idx_tgt, , drop = FALSE]
y_dis <- Y[idx_dis]
y_tgt <- Y[idx_tgt]

## prepare data for target population
Target_LA_hap1_chr <- LA_hap1_chr[idx_tgt, , drop = FALSE]
Target_LA_hap2_chr <- LA_hap2_chr[idx_tgt, , drop = FALSE]
Target_gaProp <- gaProp[idx_tgt, , drop = FALSE]
  
  
## compute PCs on training dataset
PCs_dis <-  NULL
if (use_pcs) {
  mu_dis  <- colMeans(G_dis)
  ex2_dis <- colMeans(G_dis^2)
  var_dis <- pmax(0, ex2_dis - mu_dis^2)
  sd_dis  <- sqrt(var_dis); sd_dis[sd_dis == 0] <- 1
  Z_dis <- sweep(sweep(G_dis, 2, mu_dis, `-`), 2, sd_dis, `/`)
  
  if (requireNamespace("irlba", quietly = TRUE)) {
    sv <- irlba::irlba(Z_dis, nv = n_pcs, nu = n_pcs)
    PCs_dis <- sv$u %*% diag(sv$d, nrow = n_pcs, ncol = n_pcs) # training scores = U*D
    } else {
      sv <- svd(Z_dis, nu = n_pcs, nv = n_pcs)
      PCs_dis <- sv$u %*% diag(sv$d[1:n_pcs], nrow = n_pcs, ncol = n_pcs)
      }
  colnames(PCs_dis) <- paste0("PC", seq_len(n_pcs))
  }
  

## compute single variant association analysis
## fast version
## get beta estimates

Xcov <- if (use_pcs) cbind(Intercept = 1, PCs_dis) else matrix(1, nrow(G_dis), 1, dimnames = list(NULL, "Intercept"))
qrX  <- qr(Xcov)

y_res <- qr.resid(qrX, y_dis)
G_res <- G_dis - Xcov %*% qr.coef(qrX, G_dis)

## univariate OLS on residuals: y_res ~ 0 + G_res[,j]
Gy <- as.numeric(crossprod(G_res, y_res))
G2 <- colSums(G_res^2)

beta_hat <- Gy / G2
beta_hat[!is.finite(beta_hat)] <- 0
y2 <- sum(y_res^2)
e2 <- y2 - 2 * beta_hat * Gy + (beta_hat^2) * G2
df <- nrow(G_res) - ncol(Xcov) - 1
df <- pmax(1, df)
sigma2 <- e2 / df
se <- sqrt(sigma2 / G2)
se[!is.finite(se)] <- NA_real_
z <- beta_hat / se
z[!is.finite(z)] <- 0
p <- 2 * pnorm(-abs(z))
  
## prepare final sumstat for computing PRS in validation set
sumstats <- data.table(
  rsid = val_dt$rsid,
  chrom = val_dt$chrom,
  pos = val_dt$pos,
  beta_hat = beta_hat,
  se = se,
  z = z,
  p = p,
  true_beta = val_dt$beta,
  EUR_af_alt = val_dt$EUR_af_alt,
  AFR_af_alt = val_dt$AFR_af_alt,
  EAS_af_alt = val_dt$EAS_af_alt
  )


beta_use <- sumstats$beta_hat

## only use causal variants
not_use_index <- which(sumstats$true_beta == 0)
  
beta_use[not_use_index] <- 0
beta_use[is.na(beta_use)] <- 0
PRS_tgt <- as.numeric(G_tgt %*% beta_use)

  

## compute top PCs using only target individuals
n_pcs_target_req <- 20
nv <- min(n_pcs_target_req, nrow(G_tgt), ncol(G_tgt))  # respect rank
  
## standardize target genotypes using *target* allele distribution
mu_tgt <- colMeans(G_tgt)
ex2_tgt <- colMeans(G_tgt^2)
var_tgt <- pmax(0, ex2_tgt - mu_tgt^2)
sd_tgt <- sqrt(var_tgt); sd_tgt[sd_tgt == 0] <- 1

Z_tgt_only <- sweep(sweep(G_tgt, 2, mu_tgt, `-`), 2, sd_tgt, `/`)

if (requireNamespace("irlba", quietly = TRUE)) {
  sv_t <- irlba::irlba(Z_tgt_only, nv = nv, nu = nv)
  PCs_tgt_only <- sweep(sv_t$u, 2, sv_t$d, `*`)
  svals <- sv_t$d
  } else {
    svd_t <- svd(Z_tgt_only, nu = nv, nv = nv)
    PCs_tgt_only <- sweep(svd_t$u, 2, svd_t$d[1:nv], `*`)
    svals <- svd_t$d[1:nv]
    }
  

## prepare input data for computing local ePRS
## local ancestry allele
local_anc_allele <- list(copy_1 = Target_LA_hap1_chr, copy_2 = Target_LA_hap2_chr)

merge_ancestry_alleles_to_counts <- function(local_ancestry_alleles) {
  pops <- c("EUR", "AFR", "EAS")
  n <- nrow(local_ancestry_alleles[[1]])
  p <- ncol(local_ancestry_alleles[[1]])
  counts_array <- array(0, dim = c(n, p, 3), 
                        dimnames = list(rownames(local_ancestry_alleles[[1]]),
                                        colnames(local_ancestry_alleles[[1]]),
                                        pops))
  hap1 <- local_ancestry_alleles[[1]]
  hap2 <- local_ancestry_alleles[[2]]
  for (k in seq_along(pops)) {
    pop <- pops[k]
    counts_array[,,k] <- (hap1 == pop) + (hap2 == pop)
    }
  counts_list <- lapply(seq_along(pops), function(k) {
    counts_array[,,k]
    })
  names(counts_list) <- paste0(pops, "_counts")
  return(counts_list)
  }
  
local_anc_counts <- merge_ancestry_alleles_to_counts(local_anc_allele)
  
## check final datasets for analysis
Ancestry_AF <- sumstats[, c("EUR_af_alt", "AFR_af_alt", "EAS_af_alt")]
colnames(Ancestry_AF) <- c("EUR", "AFR", "EAS")
Ancestry_AF <- as.data.frame(Ancestry_AF)

global_EPRS <- generate_global_EPRS_rPRS(PRS_tgt,
                                         Ancestry_AF,
                                         beta_use,
                                         Target_gaProp)
  
  
local_EPRS <- generate_local_EPRS_rPRS(PRS_tgt,
                                       Ancestry_AF,
                                       beta_use,
                                       local_anc_counts = local_anc_counts)
  

  
## competing methods
## ePRS
Test1 <- data.frame(Y = y_tgt, rPRS = global_EPRS$grPRS, ePRS = global_EPRS$gEPRS)

Test2 <- data.frame(Y = y_tgt, rPRS = local_EPRS$lrPRS, ePRS = local_EPRS$lEPRS)

## PRS + PCs
Test3 <- data.frame(Y = y_tgt, PRS = PRS_tgt, PCs_tgt_only[, 1:2])

Test4 <- data.frame(Y = y_tgt, PRS = PRS_tgt, PCs_tgt_only[, 1:10])

Test5 <- data.frame(Y = y_tgt, PRS = PRS_tgt, PCs_tgt_only[, 1:20])

Test6 <- data.frame(Y = y_tgt, PRS = PRS_tgt, anc1 = Target_gaProp[, 1], anc2 = Target_gaProp[, 2])

## PRS-outcome association estimates
res1 <- summary(lm(Y ~ ., data = Test1))$coef[2,"Estimate"]

res2 <- summary(lm(Y ~ ., data = Test2))$coef[2,"Estimate"]

res3 <- summary(lm(Y ~ ., data = Test3))$coef[2,"Estimate"]

res4 <- summary(lm(Y ~ ., data = Test4))$coef[2,"Estimate"]

res5 <- summary(lm(Y ~ ., data = Test5))$coef[2,"Estimate"]

res6 <- summary(lm(Y ~ ., data = Test6))$coef[2,"Estimate"]
  

c(res1, res2, res3, res4, res5, res6)



