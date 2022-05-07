library(dplyr)
library(sumstatFactors)
library(mashr)

args <- commandArgs(trailingOnly=TRUE)
fit <- readRDS(args[1])
R <- readRDS(args[2])
out_file <- args[3]


Zhat <- fit$fit$flash.fit$Y
Shat <- matrix(1, nrow = nrow(Zhat), ncol = ncol(Zhat))

data   = mash_set_data(Zhat, Shat)


m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)


Rcor <- cov2cor(R$R)
vs <- estimate_null_correlation_simple(data)
data = mash_update_data(data, V=Rcor)

U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))
U.ed = cov_ed(data, U.pca, subset=strong)

U.fit <-  cov_from_factors(t(fit$F_hat), "GFA")

U.c = cov_canonical(data)

m_ed   = mash(data, c(U.c,U.ed))
print(get_loglik(m_ed),digits = 10)

m_gfa   = mash(data, c(U.c,U.fit))
print(get_loglik(m_gfa),digits = 10)
