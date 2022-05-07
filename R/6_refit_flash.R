library(sumstatFactors)

args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]
outp <- args[2]
fit_method <- args[3]
maxiter <- as.numeric(args[4])

stopifnot(fit_method %in% c("ext", "seq"))
if(fit_method == "ext"){
    extrapolate <- TRUE
}else{
    extrapolate <- FALSE
}

fit0 <- readRDS(inp)

new_fit <- gfa_rebackfit2(fit = fit0$fit, extrapoate = extrapolate, maxiter = maxiter)
if(is.null(new_fit$fit$flash.fit$maxiter.reached)){
    fit1 <- gfa_wrapup2(new_fit$fit, nullcheck = TRUE)
}else{
    fit1 <- new_fit
}

saveRDS(fit1, outp)
