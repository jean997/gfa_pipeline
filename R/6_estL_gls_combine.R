library(dplyr)
library(purrr)
library(stringr)
library(sumstatFactors)


out <- snakemake@output[["out"]]
my_files = unlist(snakemake@input[["ests"]])

res <- map_dfr(my_files, readRDS)
saveRDS(res, file = out)
