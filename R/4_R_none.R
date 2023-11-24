

gwas_info <- snakemake@input[["gwas_info"]]
names <- readr::read_csv(gwas_info)$name
n <- length(names)
ret <- list(R = diag(n, nrow = n), names = names)
saveRDS(ret, file = snakemake@output[["out"]])
