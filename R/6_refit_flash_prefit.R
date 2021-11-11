library(sumstatFactors)

args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]
outp <- args[2]
fit_method <- args[3]
maxiter <- as.numeric(args[4])

stopifnot(fit_method %in% c("ext", "seq"))
if(fit_method == "ext"){
    fit_method <- "extrapolate"
}else{
    fit_method <- "sequential"
}

fit0 <- readRDS(inp)

new_fit <- gfa_rebackfit(fit = fit0$fit, fixed_ix = fit0$fixed_ix, method = fit_method, maxiter = maxiter)
if(is.null(new_fit$fit$flash.fit$maxiter.reached)){
    fit1 <- gfa_wrapup(new_fit$fit, new_fit$fixed_ix, nullcheck = TRUE)
}else{
    fit1 <- new_fit
}

saveRDS(fit1, outp)
