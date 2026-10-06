library(dplyr)
library(purrr)
library(stringr)
library(GFA)


out <- snakemake@output[["out"]]
R_est_file <- snakemake@input[["R"]]
params_file <- snakemake@params[["params_file"]]
max_snp <- as.numeric(snakemake@params[["max_snps"]])
seed <- snakemake@wildcards[["fs"]]
info <- readr::read_csv(snakemake@input[["gwas_info"]])

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

nms <- names(X)[grep(".z$", names(X))] %>% str_replace(".z$", "")

R <- readRDS(R_est_file)
stopifnot(all(R$names %in% nms))
z_order <- match(R$names, nms)
SS <- SS[,z_order]
Z_hat <- Z_hat[,z_order]


rownames(R$R) <- colnames(R$R) <- NULL

N <- apply(SS, 2, median)

## get case count and population prevalence
if(!"case_count" %in% names(info)){
  N_case = NULL
}else if(all(is.na(info$case_count))){
  N_case = NULL
}else{
  N_case <- info[match(R$names, info$name),]$case_count
}

if(!"pop_prev" %in% names(info)){
  pop_prev = NULL
}else if(all(is.na(info$pop_prev))){
  pop_prev = NULL
}else{
  pop_prev <- info[match(R$names, info$name),]$pop_prev
}




t <- system.time(f <- gfa_fit(Z_hat = Z_hat,
                                N = N,
                                N_case = N_case,
                                pop_prev = pop_prev,
                                R = R$R,
                                params = params))


f$snps <- snps
f$names <- R$names
f$time <- t
saveRDS(f, file=out)

