library(GFA)

inp <- snakemake@input[["inp"]]
outp <- snakemake@output[["out"]]

fit0 <- readRDS(inp)

new_fit <- gfa_rebackfit(gfa_fit = fit0,
                         params = fit0$params)
new_fit$names <- fit0$names
saveRDS(new_fit, outp)
