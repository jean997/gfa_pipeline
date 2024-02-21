library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)


# sample size affects genetic covariance and h2 but not intercept or genetic correlation
gwas_info <- read_csv(snakemake@input[["gwas_info"]])
z_files = unlist(snakemake@input[["Z"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]

ld <- purrr::map_dfr(1:22, function(c){
  read_table(ld_files[c])
})

M <- purrr:::map(1:22, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()

X <- map_dfr(z_files, function(f){
  readRDS(f) %>%
    rename(SNP = snp) %>%
    inner_join(., ld)})

names <- gwas_info$name

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()
nmsz <- str_replace(colnames(Z_hat), ".z$", "")
SS <- X %>%
  select(ends_with(".ss")) %>%
  as.matrix()
nmss <- str_replace(colnames(SS), ".ss$", "")
o <- match(nmsz, nmss)
SS <- SS[, o]
N <- apply(SS, 2, median)
if(any(is.na(N))){
  N[is.na(N)] <- gwas_info$pub_sample_size[is.na(N)]
}

R <- R_ldsc(Z_hat = Z_hat,
              ldscores = X$L2,
              ld_size = M,
              N = N,
              return_gencov = TRUE,
              return_cor = TRUE,
              make_well_conditioned = FALSE # this is done in cor_clust
              )

ret <- list(R = R$Re, Rg = R$Rg, names = nmsz)
saveRDS(ret, file=out)

