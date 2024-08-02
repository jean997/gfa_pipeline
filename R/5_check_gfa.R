library(readr)

fit <- readRDS(snakemake@input[[1]])
n <- snakemake@wildcards[["n"]]
out <- snakemake@output[["out_check"]]
sf <- snakemake@params[["success_file"]]
infile <- snakemake@input[[1]]
out_rds <- snakemake@output[["out_rds"]]


if(is.null(fit$fit$flash_fit$maxiter.reached)){
    statement <- paste0("Algorithm is converged after ", n, " rounds of backfitting.\n")
    write_lines(statement, sf)
    system(paste0("touch ", out))
}else{
    system(paste0("touch ", out))
}

