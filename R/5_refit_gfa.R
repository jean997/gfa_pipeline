library(sumstatFactors)

inp <- snakemake@input[[1]]
outp <- snakemake@output[["out"]]
params_file <- snakemake@params[["param_file"]]


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
