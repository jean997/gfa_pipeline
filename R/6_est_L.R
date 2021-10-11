library(dplyr)
library(sumstatFactors)

args <- commandArgs(trailingOnly=TRUE)
fit <- readRDS(args[1])
X <- readRDS(args[2])
R <- readRDS(args[3])
out_file <- args[4]

ntrait <- X %>%
          select(ends_with(".est")) %>%
          ncol()

nmiss <- X %>%
         select(ends_with(".est")) %>%
         is.na() %>%
         rowSums()

ix <- which(nmiss == 0)
X <- X[ix,]

nms <- names(X)[grep(".est$", names(X))]

Z_hat <- X %>%
         select(ends_with(".est")) %>%
         as.matrix()

z_order <- match(R$names, nms)
Z_hat <- Z_hat[,z_order]

R$R <- LaplacesDemon::as.symmetric.matrix(R$R)
LL <- est_L_z(Z_hat = Z_hat, R = R$R, F_hat = fit$F_hat)
nfactor <- ncol(fit$F_hat)
P <- with(LL, 2*pnorm(-abs(L_est/L_est_se)))
P <- select(X, seqnames, snp, REF, ALT, kettunen_AcAce.pos) %>% 
     rename(pos = `kettunen_AcAce.pos`) %>%
     bind_cols(., data.frame(P)) 

names(P)[-c(1:5)] <- paste0("p", 1:nfactor)
LL$P <- P

saveRDS(LL, out_file)


