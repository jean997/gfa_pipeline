library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(rlang)
library(readr)
library(purrr)
library(stringr)

source("R/format_ieu_chrom.R")

c <- as.numeric(snakemake@wildcards[["chrom"]])
gwas_info_file <- snakemake@input[["gwas_info"]]

af_thresh <- as.numeric(snakemake@params[["af_thresh"]])
sample_size_tol <- as.numeric(snakemake@params[["sample_size_tol"]])
out <- snakemake@output[["out"]]


info <- read_csv(gwas_info_file)
if(!"af" %in% names(info)){
    info$af <- NA
}

fulldat <- map(seq(nrow(info)),   function(i){
                        f <- info$raw_data_path[i]
                        if(str_ends(f, "vcf.gz") | str_ends(f, "vcf.bgz")){
                            dat <- format_ieu_chrom(f, c, af_thresh)
                        }else{
                            dat <- format_flat_chrom(f, c, af_thresh,
                                                     info$snp[i],
                                                     info$pos[i],
                                                     info$chrom[i],
                                                     info$A1[i],
                                                     info$A2[i],
                                                     info$beta_hat[i],
                                                     info$se[i],
                                                     info$p_value[i],
                                                     info$af[i],
                                                     info$sample_size[i],
                                                     as.logical(info$effect_is_or[i]))
                        }
                        if(all(is.na(dat$sample_size))){
                           dat$sample_size <- info$pub_sample_size[i]
                        }

                        if(is.finite(sample_size_tol)){
                           m <- median(dat$sample_size)
                           dat <- filter(dat, sample_size > (1-sample_size_tol)*m & sample_size < (1 + sample_size_tol)*m)
                        }
                        n <- info$name[i]
                        se_name <- as_name(paste0(n, ".se"))
                        z_name <- as_name(paste0(n, ".z"))
                        ss_name <- as_name(paste0(n, ".ss"))
                        af_name <- as_name(paste0(n, ".af"))

                        dat$sample_size[is.na(dat$sample_size)] <- as.numeric(info$pub_sample_size)
                        dat <-dat %>%  dplyr::mutate(Z = beta_hat/se) %>%
                               dplyr::rename(REF = A2, ALT = A1) %>%
                               dplyr::select(chrom, snp, REF, ALT,
                                              !!z_name := Z,
                                              !!se_name := se,
                                              !!ss_name := sample_size,
                                              !!af_name := allele_freq)
                 }) %>%
       purrr::reduce(full_join, by = c("chrom", "snp", "REF", "ALT"))

dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
if(length(dup_snps) > 0){
    fulldat <- filter(fulldat, !snp %in% dup_snps)
}



# Table of how traits are missing each SNP for LD clumping
miss <- fulldat %>%
        dplyr::select(ends_with(".z")) %>%
        is.na(.) %>%
        rowSums(.)

#nmiss <- data.frame(snp = fulldat$snp, miss = miss)
#ix <- which(miss <= nmiss_thresh)
ix <- which(miss == 0)

saveRDS(fulldat[ix,], file=out)

#saveRDS(nmiss[ix,], file=nm_out)

