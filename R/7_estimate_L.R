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


Lfit <- est_L_flash2(Z_hat = Z_hat, fit = fit$fit)
L_df <- cbind(X$snp, Lfit$L.pm, Lfit$L.psd)
names(L_df) <- c("snp", paste0(fit$names, "_pm"), paste0(fit$names, "_psd"))

saveRDS(L_df, out_file)

