library(dplyr)
library(purrr)
library(stringr)
library(MRBEE)


out <- snakemake@output[["out"]]
R <- readRDS(snakemake@input[["R"]])


pt <- as.numeric(snakemake@wildcards[["pt"]])
data_files = unlist(snakemake@input[["data"]])
max_iter = snakemake@params[["max_iter"]]
pl_pt <- as.numeric(snakemake@wildcards[["plpt"]])


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

# order input for mrbee
Rcor <- cov2cor(R$R)
ro <- c(2:p, 1)
Rcor <- Rcor[ro, ro]

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
  var_ix1 <- which(snp_info[[select_fold]]$ld_pruned == TRUE & snp_info[[select_fold]]$snp %in% var_selected)
  my_beta_hat <- beta_hat[[fit_fold]][var_ix1,]
  my_se_beta_hat <- se_beta_hat[[fit_fold]][var_ix1,]

  my_fit <- MRBEE.IMRP(by = my_beta_hat[, 1],
                     bX = my_beta_hat[, -1, drop = FALSE],
                     byse = my_se_beta_hat[, 1],
                     bXse = my_se_beta_hat[, -1, drop = FALSE],
                     Rxy = Rcor,
                     pv.thres = pl_pt,
                     FDR = F,
                     max.iter = max_iter)

  return(my_fit)
})

beta_cf <- map(fits, function(x){x$theta}) %>%
  Reduce(`+`, .)
beta_cf <- beta_cf/2

se_cf <- map(fits, function(x){diag(x$covtheta)}) %>%
  Reduce(`+`, .)
se_cf <- sqrt(se_cf)/2


ret = list(fits = fits,
           beta_cf = beta_cf,
           se_cf = se_cf,
           outcome = nms[1],
           exposures = nms[-1])


saveRDS(ret, file=out)

