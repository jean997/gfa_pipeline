library(dplyr)
library(purrr)
library(stringr)
library(esmr)


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

# select fold 2, fit on fold 1
var_selected <- dplyr::filter(snp_info[[2]], ld_pruned == TRUE) %>%
    dplyr::arrange(psel_fold) %>%
    dplyr::pull(snp)
var_selected <- var_selected[1:n]
var_ix1 <- which(snp_info[[2]]$ld_pruned == TRUE)
var_ix2 <- which(snp_info[[1]][var_ix1, "snp"] %in% var_selected)

fit1 <- esmr(beta_hat_Y = beta_hat[[1]][var_ix1,1],
            se_Y = se_beta_hat[[1]][var_ix1,1],
            beta_hat_X = beta_hat[[1]][var_ix1,2:ntrait],
            se_X = se_beta_hat[[1]][var_ix1, 2:ntrait],
            variant_ix =var_ix2,
            R = R$R,
            max_iter = max_iter)
fit1$names <- nms

# select fold 1, fit on fold 2
var_selected <- dplyr::filter(snp_info[[1]], ld_pruned == TRUE) %>%
  dplyr::arrange(psel_fold) %>%
  dplyr::pull(snp)
var_selected <- var_selected[1:n]
var_ix1 <- which(snp_info[[1]]$ld_pruned == TRUE)
var_ix2 <- which(snp_info[[2]][var_ix1, "snp"] %in% var_selected)

fit2 <- esmr(beta_hat_Y = beta_hat[[2]][var_ix1,1],
             se_Y = se_beta_hat[[2]][var_ix1,1],
             beta_hat_X = beta_hat[[2]][var_ix1,2:ntrait],
             se_X = se_beta_hat[[2]][var_ix1, 2:ntrait],
             variant_ix =var_ix2,
             R = R$R,
             max_iter = max_iter)
fit2$names <- nms

beta_cf <- (fit1$beta$beta_m + fit2$beta$beta_m)/2
se_cf <- sqrt(fit1$beta$beta_s^2 + fit2$beta$beta_s^2)/2

ret = list(fit1 = fit1,
           fit2 = fit2,
           beta_cf = beta_cf,
           se_cf = se_cf, 
           outcome = nms[1], 
           exposures = nms[-1])


saveRDS(ret, file=out)

