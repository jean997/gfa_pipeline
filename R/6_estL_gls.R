library(dplyr)
library(purrr)
library(stringr)
library(GFA)


out <- snakemake@output[["out"]]
data <- readRDS(snakemake@input[["z_file"]])
gfafit <- readRDS(snakemake@input[["gfa_file"]])
R <- readRDS(snakemake@input[["R"]])
stopifnot(all(R$names == gfafit$names))
R <- R$R

Z_hat <- data %>%
  select(ends_with(".z")) %>%
  as.matrix()

dat_names <- str_replace(colnames(Z_hat), ".z$", "")
ix <- which(dat_names %in% gfafit$names)
Z_hat <- Z_hat[,ix]
dat_names <- dat_names[ix]

o <- match(dat_names, gfafit$names)
Z_hat <- Z_hat[,o]
S <- matrix(1, nrow = nrow(Z_hat), ncol = ncol(Z_hat))
gls_sol <- GFA:::loadings_gls(Z_hat, S, R, gfafit$F_hat)
gls_zscores <- gls_sol$L/gls_sol$S
gls_pvals <- 2*pnorm(-abs(gls_zscores))
nf <- ncol(gfafit$F_hat)
res <- data.frame(cbind(gls_zscores, gls_pvals))
names(res) <- c(paste0("factor", 1:nf, ".z"), paste0("factor", 1:nf, ".p"))
res <- bind_cols(data[,c("chrom", "snp", "REF", "ALT")], res)
saveRDS(res, file = out)
