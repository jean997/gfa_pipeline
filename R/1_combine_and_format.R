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
nmiss_thresh <- as.numeric(snakemake@params[["nmiss_thresh"]])
af_thresh <- as.numeric(snakemake@params[["af_thresh"]])
out <- snakemake@output[["out"]]


info <- read_csv(gwas_info_file)
if(!af %in% names(info)) info$af <- NA

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

                        n <- info$name[i]
                        pos_name <- as_name(paste0(n, ".pos"))
                        z_name <- as_name(paste0(n, ".z"))
                        ss_name <- as_name(paste0(n, ".ss"))

                        dat$sample_size[is.na(dat$sample_size)] <- as.numeric(info$pub_sample_size)
                        dat <-dat %>%  mutate(Z = beta_hat/se) %>%
                               rename(REF = A2, ALT = A1) %>% 
                               select(chrom, snp, REF, ALT,
                                              !!pos_name := pos,
                                              !!z_name := Z,
                                              !!ss_name := sample_size)
                 }) %>%
       purrr::reduce(full_join, by = c("chrom", "snp", "REF", "ALT"))

dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
if(length(dup_snps) > 0){
    fulldat <- filter(fulldat, !snp %in% dup_snps)
}



# Save table of how traits are missing each SNP for LD clumping
miss <- fulldat %>%
        select(ends_with(".z")) %>%
        is.na(.) %>%
        rowSums(.)

#nmiss <- data.frame(snp = fulldat$snp, miss = miss)

ix <- which(miss <= nmiss_thresh)

saveRDS(fulldat[ix,], file=out)

#saveRDS(nmiss[ix,], file=nm_out)

