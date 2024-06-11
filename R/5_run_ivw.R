library(dplyr)
library(purrr)
library(esmr)

args <- commandArgs(trailingOnly=TRUE)

out <- snakemake@output[["out"]]
R_est_file <- snakemake@input[["R"]]
pt <- snakemake@wildcards[["pt"]]
pl_pt <- snakemake@wildcards[["pt"]]
z_files = unlist(snakemake@input[["Z"]])
#seed <- snakemake@wildcards[["fs"]]
#set.seed(seed)



# Read in data
X <- map_dfr(z_files, readRDS)

ntrait <- X %>%
  select(ends_with(".z")) %>%
  ncol()

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()

se_hat <- X %>%
  select(ends_with(".se")) %>%
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
#
# beta_hat <- data.frame(beta_hat)
# se_hat <- data.frame(se_hat)

pmin <- apply(pvals[,-1, drop = F], 1, min)
ix <- which(pmin < pt)

exp <- as.matrix(beta_hat[ix, 2:p])
colnames(exp) <- nms[-1]
hdat <-  list(exposure_beta = beta_hat[, 2:p],
              exposure_pval = pvals[ix, 2:p],
              exposure_se = se_hat[ix,2:p],
              outcome_beta = beta_hat[ix,1,drop = FALSE],
              outcome_pval = pvals[ix,1, drop = FALSE],
              outcome_se = se_hat[ix,1, drop = FALSE],
              expname = data.frame(id.exposure = nms, exposure = nms),
              outname = data.frame(id.outcome = nms[1], outcome = nms[1]))

t1 <- system.time(res1 <- mv_multiple(hdat,
                                      instrument_specific = FALSE))

t2 <- system.time(res2 <- mv_multiple(hdat,
                                      instrument_specific = TRUE))

res1$time <- t1
res2$time <- t2
res <- list(res1, res2)

saveRDS(res, file=out)


