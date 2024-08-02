library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)


# sample size affects genetic covariance and h2 but not intercept or genetic correlation
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

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()
nms <- str_replace(colnames(Z_hat), ".z$", "")

R <- R_ldsc(Z_hat = Z_hat,
            ldscores = X$L2,
            ld_size = M,
            N = rep(1, ncol(Z_hat)),
            return_gencov = TRUE,
            make_well_conditioned = FALSE # not needed we only need Rg
)

saveRDS(R, file=out)

