library(flashier)
library(readr)

args <- commandArgs(trailingOnly=TRUE)
fit <- readRDS(args[1])
n <- args[2]
out <- args[3]
sf <- args[4]

if(is.null(fit$fit$flash.fit$maxiter.reached)){
    statement <- paste0("Algorithm is converged after ", n, " rounds of backfitting.\n")
    write_lines(statement, sf)
    system(paste0("touch ", out))
}else{
    system(paste0("touch ", out))
}
