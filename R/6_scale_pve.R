library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(readr)
library(purrr)
library(stringr)

fit_file <- snakemake@input[["in_rds"]]
gwas_info_file <- snakemake@params[["gwas_info"]]
out <- snakemake@output[["out"]]
mode <- snakemake@wildcards[["mode"]]

info <- read_csv(gwas_info_file)
fit <- readRDS(fit_file)

if(mode == "z-score"){
    samplesize <- map(seq(nrow(info)),   function(i){
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
                 }) %>% unlist()
}else{
    samplesize <- rep(1, nrow(info))
}
ss <- data.frame(name = info$name, ss = samplesize)
s <- ss$ss[match(str_replace(fit$name, ".z", ""), ss$name)]
F_scale <- t(t(fit$F_hat)/sqrt(s))
F_scale <- apply(F_scale, 2, function(f){ f/sqrt(sum(f^2))})
fit$F_hat_effect_scale <- F_scale
fit$sample_size <- s
fit$pve <- sumstatFactors:::pve_by_trait(Lhat=fit$L_hat, Fhat=fit$F_hat, tau=fit$fit$flash_fit$tau, sample_size=s)



saveRDS(fit, file = out)




