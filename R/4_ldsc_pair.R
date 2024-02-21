library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)


t1 <- snakemake@wildcards[["name1"]]
t2 <- snakemake@wildcards[["name2"]]
l2_dir <- snakemake@params[["l2_dir"]]
gwas_info <- read_csv(snakemake@input[["gwas_info"]])
z_files = unlist(snakemake@input[["Z"]])
out <- snakemake@output[["out"]]

ld <- purrr::map_dfr(1:22, function(c){
  read_table(paste0(l2_dir, c, ".l2.ldscore.gz"))
})

M <- purrr:::map(1:22, function(c){
  read_lines(paste0(l2_dir, c, ".l2.M_5_50"))
}) %>% unlist() %>% as.numeric() %>% sum()

#bpath <- system("which bcftools", intern  = TRUE)
#bpath <- "/sw/spack/bio/pkgs/gcc-10.3.0/bcftools/1.12-g4b275ez/bin/bcftools"
#set_bcftools(bpath)

#d1 <- query_gwas(vcf = f1, chrompos = paste0(ld$CHR, ":", ld$BP)) %>% vcf_to_tibble()

X <- map_dfr(z_files, readRDS)


if(t1 == t2){
  cols <- c("snp", paste0(t1, ".z"), paste0(t1, ".ss"))
  X <- select(X, all_of(cols))
  names(X) <- c("SNP", "z1", "ss1")
  ss1 <- median(X$ss1)
  if(is.na(ss1)){
    i <- which(str_detect(t1, gwas_info$name))
    ss1 <- gwas_info$pub_sample_size[i]
  }

  full_dat <- X %>%
    inner_join(., ld)

  h2 <- snp_ldsc(ld_score = full_dat$L2,
                    ld_size = M,
                    chi2 = full_dat$z1^2,
                    sample_size = ss1,
                    blocks = NULL)

  saveRDS(h2, file= out)
}else{
  cols <- c("snp", paste0(c(t1, t2), ".z"), paste0(c(t1, t2), ".ss"))
  X <- select(X, all_of(cols))
  names(X) <- c("SNP", "z1", "z2", "ss1", "ss2")
  ss1 <- median(X$ss1)
  ss2 <- median(X$ss2)
  if(is.na(ss1)){
    i <- which(str_detect(t1, gwas_info$name))
    ss1 <- gwas_info$pub_sample_size[i]
  }
  if(is.na(ss2)){
        i <- which(str_detect(t2, gwas_info$name))
        ss2 <- gwas_info$pub_sample_size[i]
  }

  full_dat <- inner_join(X, ld)

  rg <- ldsc_rg(ld_score = full_dat$L2,
                  ld_size = M,
                  sample_size_1 = ss1,
                  sample_size_2 = ss2,
                  z1 = full_dat$z1,
                  z2 = full_dat$z2,
                  blocks = NULL, h2_se = FALSE)

  saveRDS(rg, file= out)
}
