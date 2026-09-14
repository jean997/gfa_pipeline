library(dplyr)
library(purrr)
library(stringr)
library(TwoSampleMR)


out <- snakemake@output[["out"]]
R <- readRDS(snakemake@input[["R"]])

pt <- as.numeric(snakemake@wildcards[["pt"]])
data_files = unlist(snakemake@input[["data"]])
max_iter = snakemake@params[["max_iter"]]
#seed <- snakemake@wildcards[["fs"]]
#set.seed(seed)



# Read in data
split_data <- lapply(data_files, readRDS)
#nfolds <- length(split_data) # should be 2

ntrait <- ncol(split_data[[1]][[1]]$beta_hat)
nms <- split_data[[1]][[1]]$names

snp_info <- purrr::map(1:2, function(i){
  map(split_data, i) %>% map_dfr("snp_info")
})

beta_hat <- purrr::map(1:2, function(i){
  map(split_data, i) %>% map_dfr(function(x){data.frame(x[["beta_hat"]])})
})

se_beta_hat <- purrr::map(1:2, function(i){
  map(split_data, i) %>% map_dfr(function(x){data.frame(x[["se_beta_hat"]])})
})

n <- sum(snp_info[[1]]$psel_full < pt & snp_info[[1]]$full_pruned == TRUE) # should be the same for fold 1 or 2
p <- ncol(beta_hat[[1]])

# select fold 2, fit on fold 1
fits <- lapply(1:2, function(i){
  fit_fold <- i
  if(fit_fold == 1){
    select_fold <- 2
  }else if(fit_fold == 2){
    select_fold <- 1
  }

  var_selected <- dplyr::filter(snp_info[[select_fold]], ld_pruned == TRUE) %>%
    dplyr::arrange(psel_fold) %>%
    dplyr::pull(snp)
  var_selected <- var_selected[1:n]
  effective_pt <- filter(snp_info[[select_fold]], snp %in% var_selected) %>% pull(psel_fold) %>% max()

  var_ix1 <- which(snp_info[[select_fold]]$ld_pruned == TRUE & snp_info[[select_fold]]$snp %in% var_selected)
  my_beta_hat <- as.matrix(beta_hat[[fit_fold]][var_ix1,])
  my_se_beta_hat <- as.matrix(se_beta_hat[[fit_fold]][var_ix1,])

  Z_sel <- as.matrix(beta_hat[[select_fold]]/se_beta_hat[[select_fold]])[var_ix1,]
  pval_sel <- 2*pnorm(-abs(Z_sel))

  colnames(my_beta_hat) <- nms
  hdat <-  list(exposure_beta = my_beta_hat[, 2:p],
                 exposure_pval = pval_sel[, 2:p],
                 exposure_se = my_se_beta_hat[,2:p],
                 outcome_beta = my_beta_hat[,1],
                 outcome_pval = pval_sel[,1],
                 outcome_se = my_se_beta_hat[,1],
                 expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                 outname = data.frame(id.outcome = nms[1], outcome = nms[1]))

  fit_nspec <- mv_multiple(hdat, instrument_specific = FALSE, pval_threshold = effective_pt)
  fit_spec <- mv_multiple(hdat, instrument_specific = TRUE, pval_threshold = effective_pt)

  return(list("instrument_specific" = fit_spec, "not_instrument_specific"= fit_nspec))
})

beta_cf_spec <- map(fits, function(x){x$instrument_specific$result$b}) %>%
  Reduce(`+`, .)
beta_cf_spec <- beta_cf_spec/2

se_cf_spec <- map(fits, function(x){x$instrument_specific$result$se^2}) %>%
  Reduce(`+`, .)
se_cf_spec <- sqrt(se_cf_spec)/2


beta_cf_nspec <- map(fits, function(x){x$not_instrument_specific$result$b}) %>%
  Reduce(`+`, .)
beta_cf_nspec <- beta_cf_nspec/2

se_cf_nspec <- map(fits, function(x){x$not_instrument_specific$result$se^2}) %>%
  Reduce(`+`, .)
se_cf_nspec <- sqrt(se_cf_nspec)/2

ret = list(fits = fits,
           beta_cf_spec = beta_cf_spec,
           se_cf_spec = se_cf_spec,
           beta_cf_nspec = beta_cf_nspec,
           se_cf_nspec = se_cf_nspec,
           outcome = nms[1],
           exposures = nms[-1])


saveRDS(ret, file=out)

