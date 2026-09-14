library(dplyr)
library(purrr)
library(stringr)
library(esmr)


out <- snakemake@output[["out"]]
R <- readRDS(snakemake@input[["R"]])
Rcor <- cov2cor(R$R)

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
  selection_pvalue <- filter(snp_info[[select_fold]], snp %in% var_selected) %>% pull(psel_fold)
  effective_pt <- max(selection_pvalue)

  var_ix1 <- which(snp_info[[select_fold]]$ld_pruned == TRUE & snp_info[[select_fold]]$snp %in% var_selected)
  if(length(var_ix1) > 5000){
    var_ix1 <- sample(var_ix1, size = 5000, replace = FALSE)
  }


  my_beta_hat <- beta_hat[[fit_fold]][var_ix1,]
  my_se_beta_hat <- se_beta_hat[[fit_fold]][var_ix1,]

  grapple_dat <- data.frame(cbind(my_beta_hat, my_se_beta_hat))
  names(grapple_dat) <- c("gamma_out", paste0("gamma_exp", 1:(p-1)),
                          "se_out", paste0("se_exp", 1:(p-1)))

  grapple_dat$selection_pvals <- selection_pvalue


  my_fit <- grappleRobustEst(data = grapple_dat,
                            plot.it =FALSE,
                            p.thres = effective_pt,
                            cor.mat = Rcor,
                            niter = max_iter)

  return(my_fit)
})

beta_cf <- map(fits, function(x){x$beta.hat}) %>%
  Reduce(`+`, .)
beta_cf <- beta_cf/2

se_cf <- map(fits, function(x){diag(x$beta.var)}) %>%
  Reduce(`+`, .)
se_cf <- sqrt(se_cf)/2


ret = list(fits = fits,
           beta_cf = beta_cf,
           se_cf = se_cf,
           outcome = nms[1],
           exposures = nms[-1])


saveRDS(ret, file=out)

