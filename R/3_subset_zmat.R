library(dplyr)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)
normbeta <- readRDS(args[1])
keep_snps <- readRDS(args[2])
out_data <- args[3]


nms <- names(normbeta)[grep(".est$", names(normbeta))]
n <- length(nms)

X <- normbeta %>%
      filter(snp %in% keep_snps$snp) 

saveRDS(X, file=out_data)

