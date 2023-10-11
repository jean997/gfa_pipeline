library(dplyr)
library(purrr)
library(stringr)
library(sumstatFactors)


out <- snakemake@output[["out"]]
mode <- snakemake@wildcards[["mode"]]
R_est_file <- snakemake@input[["R"]]
params_file <- snakemake@params[["params_file"]]
max_snp <- snakemake@params[["max_snps"]]
seed <- snakemake@wildcards[["fs"]]

z_files = unlist(snakemake@input[["Z"]])



set.seed(seed)



if(params_file == "default"){
  params <- gfa_default_parameters()
}else{
  params <- readRDS(params_file)
}


# Read in data
X <- map_dfr(z_files, readRDS)

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

if(str_ends(R_est_file, "none_R.txt")){
  R <- list(names = nms, R = diag(length(nms)))
}else{
  R <- readRDS(R_est_file)
  stopifnot(all(R$names %in% nms))
  z_order <- match(R$names, nms)
  SS <- SS[,z_order]
  Z_hat <- Z_hat[,z_order]
}




if(mode == "z-score"){
  f <- gfa_fit(Z_hat = Z_hat, R = R$R, params = params)
  N <- apply(SS, 2, median)
  f$F_hat_scaled <- t(t(f$F_hat)/sqrt(N))
}else{
  N <- apply(SS, 2, median)
  B_std <- t( t(Z_hat)/sqrt(N))
  f <- gfa_fit(B_std = B_std, N = N, R = R$R, params = params)
}

f$snps <- snps
f$names <- R$names

saveRDS(f, file=out)

