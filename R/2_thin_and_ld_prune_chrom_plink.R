library(dplyr)
library(ieugwasr)
library(readr)
library(stringr)
library(purrr)

source("R/thin_gwas.R")
X <- readRDS(snakemake@input[["zmat"]])
R <- readRDS(snakemake@input[["R"]])
region_bed <- read_tsv(snakemake@input[["region_file"]]) # should be a bed file
X_pruned <- readRDS(snakemake@input[["full_pruned"]])


r2_thresh <- as.numeric(snakemake@wildcards[["r2_thresh"]])
clump_kb <- snakemake@wildcards[["kb"]]
chrom <- as.numeric(snakemake@wildcards[["chrom"]])

ref_path  <- snakemake@params[["ref_path"]]
bim <- read_tsv(paste0(ref_path, ".bim"), col_names = c("chr", "snp", "cm", "pos", "A1", "A2"))

#out_full <- snakemake@output[["out_full"]]
out_split <- snakemake@output[["out_split"]]

pthresh <- as.numeric(snakemake@params[["pthresh"]])
pthresh_thin <- as.numeric(snakemake@params[["pthresh_thin"]])
thin_eps <- as.numeric(snakemake@params[["thin_eps1"]])

is_mvmr <- as.numeric(snakemake@params[["is_mvmr"]])

#if(!pthresh_thin >= pthresh){
#  stop("pthresh_thin must be greater than or equal to pthresh.\n")
#}


Z_hat <- X %>%
  dplyr::select(ends_with(".z")) %>%
  as.matrix()

pos <- X %>%
  dplyr::select(ends_with(".pos")) %>%
  as.matrix()
pos <- pos[,1] # assuming all the gwas are on the same position system, should check this is true eventually

se_beta_hat_orig <- X %>%
  dplyr::select(ends_with(".se")) %>%
  as.matrix()

beta_hat_orig <- Z_hat*se_beta_hat_orig

nms <- names(X)[grep(".z$", names(X))] %>% str_replace(".z$", "")
nms2 <- names(X)[grep(".se$", names(X))] %>% str_replace(".se$", "")
stopifnot(all(nms == nms2))
stopifnot(all(R$names  == nms))

#####################################################
## Reduce to snps in LD panel passing pthresh_thin ##
#####################################################

if(is_mvmr == 1){
  zmax <- apply(abs(Z_hat[, -1, drop = FALSE]), 1, function(x){max(x, na.rm=T)})
  psel_full <- 2*pnorm(-abs(zmax))
}else{
  zmax <- apply(abs(Z_hat), 1, function(x){max(x, na.rm=T)})
  psel_full <- 2*pnorm(-abs(zmax))
}
# approximation to reduce computational burden, pthresh_thin should be substantially bigger than the selection threshold.
ix <- which(psel_full < pthresh_thin & X$snp %in% bim$snp)
Z_hat <- Z_hat[ix,]
beta_hat_orig <- beta_hat_orig[ix,]
se_beta_hat_orig <- se_beta_hat_orig[ix,]
pos <- pos[ix]
snp <- X$snp[ix]
psel_full <- psel_full[ix]

#############################
## LD Pruning on Full data ##
#############################

# clump_res_full <- ld_clump(dat = data.frame(rsid = snp,
#                                             pval = psel_full),
#                            clump_r2 = r2_thresh,
#                            clump_p = pthresh,
#                            clump_kb = clump_kb,
#                            plink_bin = genetics.binaRies::get_plink_binary(),
#                            bfile = ref_path)
#
# X_pruned <- X[X$snp %in% clump_res_full$rsid,]
# saveRDS(X_pruned, file=out_full)
# full_pruned <- snp %in% clump_res_full$rsid

###########################
## Perform data thinning ##
###########################

# break snps into LD blocks
r <- dplyr::filter(region_bed, chr == paste0("chr", chrom))
snp_region <- findInterval(pos, r$start)
unique_regions <- unique(snp_region)

R_LD <- lapply(unique_regions, function(i){
  cat("Region ", i, "of ", max(unique_regions), "\n")
  snps <- snp[snp_region == i]
  ld_mat <- ieugwasr::ld_matrix(variants = snps,
                                bfile = ref_path,
                                plink_bin = genetics.binaRies::get_plink_binary(),
                                with_alleles = FALSE)
  ld_mat_snps <- rownames(ld_mat)
  if(!all(ld_mat_snps == snps)){
      stopifnot(length(ld_mat_snps) == length(snps))
      stopifnot(all(ld_mat_snps %in% snps))
      o <- match(snps, ld_mat_snps)
      ld_mat <- ld_mat[o, o]
      ld_mat_snps <- rownames(ld_mat)
      stopifnot(all(ld_mat_snps == snps))
  }
  return(list("ld_mat" = ld_mat, "ld_mat_snps" = ld_mat_snps))
})
snps_in_ld_ref <- purrr::map(R_LD, "ld_mat_snps") %>% unlist()
stopifnot(all(snps_in_ld_ref == snp)) # should match exactly and be in the right order

R_LD <- purrr::map(R_LD, "ld_mat") %>% purrr::map(function(x){as.matrix(Matrix::nearPD(x)$mat)})

#beta_hat_orig <- beta_hat
#se_beta_hat_orig <- se_hat


split_data_chrom <- thin_gwas(beta_hat = beta_hat_orig,
                              se_beta_hat = se_beta_hat_orig,
                              R = R$R,
                              R_LD = R_LD,
                              eps = c(thin_eps, 1-thin_eps))


for(i in seq(length(split_data_chrom))){
  my_z_hat <- with(split_data_chrom[[i]], beta_hat/se_beta_hat)
  if(is_mvmr == 1){
    zmax <- apply(abs(my_z_hat[, -1, drop = FALSE]), 1, function(x){max(x, na.rm=T)})
    psel_fold <- 2*pnorm(-abs(zmax))
  }else{
    zmax <- apply(abs(my_z_hat), 1, function(x){max(x, na.rm=T)})
    psel_fold <- 2*pnorm(-abs(zmax))
  }

  clump_res <- ieugwasr::ld_clump(dat = data.frame(rsid = snp,
                                                   pval = psel_fold),
                                  bfile = ref_path,
                                  clump_r2 = r2_thresh,
                                  clump_p = pthresh,
                                  clump_kb = clump_kb,
                                  plink_bin = genetics.binaRies::get_plink_binary())
  split_data_chrom[[i]]$snp_info <- data.frame(snp = snp,
                                               psel_fold = psel_fold,
                                               psel_full = psel_full,
                                               ld_pruned = snp %in% clump_res$rsid,
                                               full_pruned = snp %in% X_pruned$snp)
  split_data_chrom[[i]]$names <- nms
}

saveRDS(split_data_chrom, file = out_split)





