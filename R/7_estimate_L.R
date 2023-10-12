library(dplyr)
library(sumstatFactors)

args <- commandArgs(trailingOnly=TRUE)
fit <- readRDS(args[1])
X <- readRDS(args[2])
out_file <- args[3]

ntrait <- X %>%
  select(ends_with(".z")) %>%
  ncol()

nmiss <- X %>%
  select(ends_with(".z")) %>%
  is.na() %>%
  rowSums()

ix <- which(nmiss == 0)
X <- X[ix,]

nms <- names(X)[grep(".z$", names(X))]

Z_hat <- X %>%
         select(ends_with(".z")) %>%
         as.matrix()


z_order <- match(fit$names, nms)
Z_hat <- Z_hat[,z_order]

nfct <- ncol(fit$F_hat)
Lfit <- est_L_flash2(Z_hat = Z_hat, fit = fit$fit)
L_df <- cbind( Lfit$L.pm[,-fit$fixed_ix], Lfit$L.psd[, -fit$fixed_ix],
               Lfit$L.lfsr[, -fit$fixed_ix]) %>% data.frame()
names(L_df) <- c(paste0("fct_", 1:nfct, "_pm"),
                 paste0("fct_", 1:nfct, "_psd"),
                 paste0("fct_", 1:nfct, "_lfsr"))
L_df$snp <- X$snp
L_df$chr <- X$chr

saveRDS(L_df, out_file)

