library(dplyr)
library(purrr)
library(stringr)
library(GFA)



p_thresh <- as.numeric(snakemake@wildcards[["pt"]])
out <- snakemake@output[["out"]]
z_files = unlist(snakemake@input[["Z"]])
cond_num = as.numeric(snakemake@params[["cond_num"]])

if(!is.finite(cond_num)){
  mwc <- FALSE
}else{
  mwc <- TRUE
}

# Read in data
X <- map_dfr(z_files, readRDS)

ntrait <- X %>%
  select(ends_with(".z")) %>%
  ncol()

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()

nms <- colnames(Z_hat) %>% stringr::str_replace(".z$", "")
Rpt <- R_pt(B_hat = Z_hat,
            S_hat = matrix(1, nrow = nrow(Z_hat), ncol = ncol(Z_hat)),
            p_val_thresh = p_thresh,
            return_cor = TRUE,
            make_well_conditioned = mwc,
            cond_num = cond_num
            )

ret <- list(R = Rpt, names = nms)
saveRDS(ret, file=out)
