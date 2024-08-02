library(dplyr)
library(purrr)
library(MRBEE)
library(stringr)

out <- snakemake@output[["out"]]
R_est_file <- snakemake@input[["R"]]
pt <- as.numeric(snakemake@wildcards[["pt"]])
pl_pt <- as.numeric(snakemake@wildcards[["plpt"]])
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
pvals <- 2*pnorm(-abs(Z_hat))

nms <- names(X)[grep(".z$", names(X))] %>% str_replace(".z$", "")
nms2 <- names(X)[grep(".se$", names(X))] %>% str_replace(".se$", "")
stopifnot(all(nms == nms2))


R <- readRDS(R_est_file)
stopifnot(all(R$names  == nms))
# z_order <- match(R$names, nms)
# se_hat <- se_hat[,z_order]
# beta_hat <- beta_hat[,z_order]
p <- ncol(beta_hat)

Rcor <- cov2cor(R$R)


ro <- c(2:p, 1)
Rcor <- Rcor[ro, ro]

pmin <- apply(pvals[,-1, drop  = F], 1, min)
ix <- which(pmin < pt)

t1 <- system.time(fit <- MRBEE.IMRP(by = beta_hat[ix, 1],
                                    bX = beta_hat[ix, -1, drop = FALSE],
                                    byse = se_hat[ix, 1],
                                    bXse = se_hat[ix, -1, drop = FALSE],
                                    Rxy = Rcor,
                                    pv.thres = pl_pt,
                                    FDR = F,
                                    max.iter = max_iter))


fit$time <- t1
fit$names <- nms
fit$z <- fit$theta/sqrt(diag(fit$covtheta))
fit$p <- 2*pnorm(-abs(fit$z))
cat("four\n")
saveRDS(fit, file=out)

