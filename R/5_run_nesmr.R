library(dplyr)
library(purrr)
library(stringr)
library(esmr)


out <- snakemake@output[["out"]]
R_est_file <- snakemake@input[["R"]]
pt <- as.numeric(snakemake@wildcards[["pt"]])
z_files = unlist(snakemake@input[["Z"]])
max_iter = snakemake@params[["max_iter"]]
template_file = snakemake@input[["template"]]
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

Rcor <- cov2cor(R$R)
p <- ncol(beta_hat)

template <- readRDS(template_file)
stopifnot(nrow(template) == nrow(Rcor))

pmin <- apply(pvals, 1, min)
ix <- which(pmin < pt)

t <- system.time(
  fit <- esmr(beta_hat_X = beta_hat,
              se_X = se_hat,
              R = Rcor,
              variant_ix = ix,
              direct_effect_template = template,
              G = diag(p),
              max_iter = max_iter))
fit$time <- t
fit$names <- nms
fit$likelihood <- esmr:::log_py(fit)

saveRDS(fit, file=out)

