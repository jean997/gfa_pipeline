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

saveRDS(new_fit, outp)
