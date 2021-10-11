library(dplyr)
library(purrr)
library(readr)

args <- commandArgs(trailingOnly=TRUE)
X <- readRDS(args[1])
l2 <- read_table(args[2])
out_summ <- args[3]


nms <- names(X)[grep(".z$", names(X))]
n <- length(nms)

X <- select(l2, SNP, L2) %>%
     rename(snp  = SNP) %>%
     left_join(., X, by = "snp")


res <- expand.grid(n1 = seq(n), n2 = seq(n)) %>%
            filter(n1 <= n2) 
res$name1 <- nms[res$n1]
res$name2 <- nms[res$n2]

summ <- map2_df(res$name1, res$name2, function(n1, n2){
             y <- X[[n1]]*X[[n2]]
             x <- X$L2
             ix <- which(!is.na(x) & !is.na(y))
             x <- x[ix]
             y <- y[ix]
             xty <- sum(x*y)
             m <- length(ix)
             xtx <- sum(x^2)
             xsum <- sum(x)
             ysum <- sum(y)
             return(list(xty = xty, xtx = xtx, xsum = xsum, ysum = ysum, m  = m, n1 = n1, n2  = n2))
             })

saveRDS(summ, file=out_summ)




