

gwas_info <- snakemake@input[["gwas_info"]]
names <- read_csv(gwas_info)$name
n <- length(names)
R <- list(R = diag(n, nrow = n), names = names)
