library(dplyr)
library(purrr)
library(readr)


out <- snakemake@output[["out"]]
gwas_info <- snakemake@input[["gwas_info"]]
root <- snakemake@params[["root"]]


names <- read_csv(gwas_info)$name

df <- expand.grid(n1 = names, n2 = names) %>%
      mutate(file1 = paste0(root, "ldsc.", n1, "__", n2, ".RDS"),
             file2 = paste0(root, "ldsc.", n2, "__", n1, ".RDS"))

df$exists1 <- file.exists(df$file1)
df$exists2 <- file.exists(df$file2)

df <- df %>% mutate(file = case_when(exists1 ~ file1, exists2 ~ file2, TRUE ~ NA_character_))

df$cov <- sapply(df$file, function(f){
  x <- readRDS(f)
  x[["int"]]
})

cov_mat <- df %>%
  select(n1, n2, cov) %>%
  reshape2::dcast(n1 ~ n2)

nms <- paste0(as.vector(cov_mat$n1), ".z")
R <- as.matrix(cov_mat[,-1])
R <- cov2cor(R)

vals <- eigen(R, only.values = TRUE)
if(any(vals) < 0){
  R <- Matrix::nearPD(R, corr = TRUE, posd.tol = 1e-3)$mat |> as.matrix()
}

#eS <- eigen(R)

ret <- list(R = R, names = nms) #, eS = eS)

saveRDS(ret, file=out)
