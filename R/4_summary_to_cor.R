library(dplyr)
library(purrr)

out = snakemake@output[["out"]]
files = unlist(snakemake@input)


x <- map(files, readRDS)


df <- map(x, function(y){mutate(y, prod = unlist(prod),
                                   ss = unlist(n),
                                   n1s = unlist(n1s),
                                   n2s  = unlist(n2s))}) %>%
        reduce(full_join, by=c("n1", "n2")) %>%
        mutate(prod = rowSums(select(., starts_with("prod"))),
               n = rowSums(select(., starts_with("ss"))),
               n1s = rowSums(select(., starts_with("n1s"))),
               n2s = rowSums(select(., starts_with("n2s")))) %>%
        select(n1, n2, prod, n, n1s, n2s) %>%
        mutate(cov = 1/(n-1)*(prod - n1s*n2s/n))

# we want a symmetric matrix so we need to replicate some rows of df
df_copy <- filter(df, n1 != n2) %>%
           rename(n1c = n2, n2c = n1) %>%
           rename(n1 = n1c, n2 = n2c)

cov_mat <- bind_rows(df, df_copy)  %>%
           select(n1, n2, cov) %>%
           reshape2::dcast(n1 ~ n2)

nms <- cov_mat$n1 %>% stringr::str_replace(".z$", "")
R <- as.matrix(cov_mat[,-1])
## no cov2cor needed, should already be correlation matrix

vals <- eigen(R, only.values = TRUE)$values
if(any(vals < 0)){
  R <- Matrix::nearPD(R, corr = TRUE, posd.tol = 1e-3)$mat |> as.matrix()
}

#eS <- eigen(R)

ret <- list(R = R, names = nms) # eS = eS)

saveRDS(ret, file=out)


