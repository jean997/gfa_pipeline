library(dplyr)
library(purrr)
library(stringr)
library(sumstatFactors)

args <- commandArgs(trailingOnly=TRUE)

out <- args[1]
mode <- args[2]
R_est_file <- args[3]
params_file <- args[4]
max_snp <- args[5]
seed <- as.numeric(args[5])

nb_files = args[-c(1:5)]

set.seed(seed)



if(params_file == "default"){
  params <- gfa_default_parameters()
}else{
  params <- readRDS(params_file)
}


# Read in data
X <- map_dfr(nb_files, readRDS)

ntrait <- X %>%
          select(ends_with(".z")) %>%
          ncol()

if(nrow(X) > max_snp){
    ix <- sample(seq(nrow(X)), size = max_snp, replace = FALSE)
    X <- X[ix,]
}


Z_hat <- X %>%
         select(ends_with(".z")) %>%
         as.matrix()

SS <- X %>%
      select(ends_with(".ss")) %>%
      as.matrix()

snps <- X$snp

nms <- names(X)[grep(".z$", names(X))]

if(R_est_file == "none"){
  R <- list(names = nms, R = diag(length(nms)))
}else{
  R <- readRDS(R_est_file)
  stopifnot(all(R$names %in% nms))
  z_order <- match(R$names, nms)
  SS <- SS[,z_order]
  Z_hat <- Z_hat[,z_order]
}




if(mode == "z-hat"){
  f <- gfa_fit(Z_hat = Z_hat, R = R$R, params = params)
}else{
  N <- apply(SS, 2, median)
  B_std <- t( t(Z_hat)*N)
  f <- gfa_fit(B_std = B_std, N = N, R = R$R, params = params)
}

f$snps <- snps
f$names <- R$names

saveRDS(f, file=out)

