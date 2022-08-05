library(dplyr)
library(purrr)

library(VariantAnnotation)
library(gwasvcf)
library(readr)
library(bigsnpr)
library(stringr)


args <- commandArgs(trailingOnly=TRUE)
f1 <- args[1]
f2 <- args[2]
l2_dir <- args[3]
gwas_info <- read_csv(args[4])
out <- args[5]


ld <- purrr::map_dfr(1:22, function(c){
  read_table(paste0(l2_dir, c, ".l2.ldscore.gz"))
})

M <- purrr:::map(1:22, function(c){
  read_lines(paste0(l2_dir, c, ".l2.M_5_50"))
}) %>% unlist() %>% as.numeric() %>% sum()


set_bcftools()

d1 <- query_gwas(vcf = f1, chrompos = paste0(ld$CHR, ":", ld$BP)) %>% vcf_to_tibble()
ss1 <- median(d1$SS)
if(is.na(ss1)){
    i <- which(str_detect(f1, gwas_info$name))
    ss1 <- gwas_info$pub_sample_size[i]
}
dat1 <- d1 %>%
        mutate(z1 = ES/SE) %>%
        dplyr::select(rsid, z1)


if(f1 == f2){
  full_dat <- dat1 %>%
    dplyr::rename(SNP = rsid) %>%
    inner_join(., ld)

  h2 <- snp_ldsc(ld_score = full_dat$L2,
                    ld_size = M,
                    chi2 = full_dat$z1^2,
                    sample_size = ss1,
                    blocks = NULL)

  saveRDS(h2, file= out)
}else{
  d2 <- query_gwas(vcf = f2, chrompos = paste0(ld$CHR, ":", ld$BP)) %>% vcf_to_tibble()
  ss2 <- median(d2$SS)
  dat2 <- d2 %>%
        mutate(z2 = ES/SE) %>%
        dplyr::select(rsid, z2)

   if(is.na(ss2)){
        i <- which(str_detect(f2, gwas_info$name))
        ss2 <- gwas_info$pub_sample_size[i]
   }

  full_dat <- inner_join(dat1, dat2) %>%
            dplyr::rename(SNP = rsid) %>%
            inner_join(., ld)

  rg <- snp_ldsc_rg(ld_score = full_dat$L2,
                  ld_size = M,
                  sample_size_1 = ss1,
                  sample_size_2 = ss2,
                  z1 = full_dat$z1,
                  z2 = full_dat$z2,
                  blocks = NULL, h2_se = FALSE)

  saveRDS(rg, file= out)
}
