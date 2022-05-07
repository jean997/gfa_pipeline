library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(readr)
library(purrr)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
fit_file <- args[1]
gwas_info_file <- args[2]
dir <- args[3]
out <- args[4]

info <- read_csv(gwas_info_file)
fit <- readRDS(fit_file)

samplesize <- map(seq(nrow(info)),   function(i){
                        #f <- paste0(dir, info$name[i], ".vcf.bgz")
                        #n <- info$name[i]
                        ss <- info$pub_sample_size[i]
                        if(is.na(ss)) ss <- info$pub_cases[i] + info$pub_controls[i]
                        return(ss)
                        dat <- map_dfr(1:22, function(c){
                            cat(c, " ")
                            v <- query_chrompos_file(paste0(c, ":1-536870911"), f)
                            d <- vcf_to_tibble(v) %>%
                                   dplyr::rename(snp=rsid) %>%
                                   dplyr::filter(snp %in% fit$snps)
                             return(d)

                        })
                        dat$SS[is.na(dat$SS)] <- ss
                        sample_size <- median(dat$SS)

                        return(sample_size)
                 })
ss <- data.frame(name = info$name, ss = unlist(samplesize))
s <- ss$ss[match(str_replace(fit$name, ".z", ""), ss$name)]
F_scale <- fit$F_hat/s
F_scale <- apply(F_scale, 2, function(f){ f/sqrt(sum(f^2))})
fit$F_scale <- F_scale
saveRDS(fit, file = out)




