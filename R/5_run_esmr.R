library(dplyr)
library(purrr)
library(stringr)
library(esmr)


out <- snakemake@output[["out"]]
R_est_file <- snakemake@input[["R"]]
pt <- as.numeric(snakemake@wildcards[["pt"]])
z_files = unlist(snakemake@input[["Z"]])
max_iter = snakemake@params[["max_iter"]]
#seed <- snakemake@wildcards[["fs"]]
#set.seed(seed)



# Read in data
X <- map_dfr(z_files, readRDS)

ntrait <- X %>%
  dplyr::select(ends_with(".z")) %>%
  ncol()

Z_hat <- X %>%
  dplyr::select(ends_with(".z")) %>%
  as.matrix()

se_hat <- X %>%
  dplyr::select(ends_with(".se")) %>%
  as.matrix()

beta_hat <- Z_hat*se_hat

nms <- names(X)[grep(".z$", names(X))] %>% str_replace(".z$", "")
nms2 <- names(X)[grep(".se$", names(X))] %>% str_replace(".se$", "")
stopifnot(all(nms == nms2))


R <- readRDS(R_est_file)
stopifnot(all(R$names  == nms))
# z_order <- match(R$names, nms)
# se_hat <- se_hat[,z_order]
# beta_hat <- beta_hat[,z_order]

Rcor <- cov2cor(R$R)
p <- ncol(beta_hat)

t <- system.time(
  fit <- esmr(beta_hat_Y = beta_hat[,1],
              se_Y = se_hat[,1],
              beta_hat_X = beta_hat[,2:p],
              se_X = se_hat[, 2:p],
              R = Rcor,
              pval_thresh = pt,
              max_iter = max_iter))
fit$time <- t
fit$names <- nms
saveRDS(fit, file=out)

