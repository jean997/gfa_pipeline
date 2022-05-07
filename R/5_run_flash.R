library(dplyr)
library(purrr)
library(stringr)
library(sumstatFactors)

args <- commandArgs(trailingOnly=TRUE)

out <- args[1]
R_est_file <- args[2]
type <- args[3]
kmax <- as.numeric(args[4])
max_snp <- as.numeric(args[5])


fit_method <- args[6]
maxiter <- as.numeric(args[7])
seed <- as.numeric(args[8])

nb_files = args[-c(1:8)]

set.seed(seed)


stopifnot(fit_method %in% c("ext", "seq"))
if(fit_method == "ext"){
    extrapolate <- TRUE
}else{
    extrapolate <- FALSE
}



stopifnot(type =="plain" | str_starts(type, "ff-"))
if(str_starts(type, "ff")){
    max_ev <- as.numeric(str_replace(type, "ff-", ""))
    type <- "ff"
}

# Read in data
X <- map_dfr(nb_files, readRDS)

ntrait <- X %>%
          select(ends_with(".z")) %>%
          ncol()
nmiss <- X %>%
         select(ends_with(".z")) %>%
         is.na() %>%
         rowSums()


ix <- which(nmiss == 0)
if(length(ix) > max_snp){
    ix <- sample(ix, size = max_snp, replace = FALSE)
}
X <- X[ix,]


Z_hat <- X %>%
         select(ends_with(".z")) %>%
         as.matrix()

SS <- X %>%
      select(ends_with(".ss")) %>%
      as.matrix()

snps <- X$snp

nms <- names(X)[grep(".z$", names(X))]

R <- readRDS(R_est_file)
stopifnot(all(R$names %in% nms))

z_order <- match(R$names, nms)
SS <- SS[,z_order]
Z_hat <- Z_hat[,z_order]
new_names <- R$names


if(type == "plain"){
        f <- fit_ff_update(Z_hat = Z_hat, kmax = kmax,
                           extrapolate = extrapolate, max_final_iter = maxiter,
                           max_lr_percent = max_ev)
}else if(type == "ff"){
        f <- fit_ff_update(Z_hat = Z_hat, R = R$R, kmax = kmax,
                           extrapolate = extrapolate, max_final_iter = maxiter,
                           max_lr_percent = max_ev)
}

f$snps <- snps
f$names <- new_names

saveRDS(f, file=out)

