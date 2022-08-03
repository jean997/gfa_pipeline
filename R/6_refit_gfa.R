library(sumstatFactors)

args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]
outp <- args[2]
params_file <- args[3]

param_default <- gfa_default_parameters()
if(params_file == "default"){
  params <- param_default
}else{
  params <- readRDS(params_file)
  if(is.null(params$extrapolate)) params$extrapolate <- param_default$extrapolate
  if(is.null(params$max_iter)) params$max_iter <- param_default$max_iter
}

fit0 <- readRDS(inp)

new_fit <- gfa_rebackfit(fit = fit0$fit,
                         params = fit0$params)
if(is.null(new_fit$fit$flash.fit$maxiter.reached)){
    fit1 <- gfa_wrapup(new_fit$fit, nullcheck = TRUE)
    fit1$params <- fit0$params
}else{
    fit1 <- new_fit
}

saveRDS(fit1, outp)
