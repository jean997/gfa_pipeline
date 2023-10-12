library(dplyr)
library(purrr)

X <- readRDS(snakemake@input[["zmat"]])
p_thresh <- as.numeric(snakemake@wildcards[["pt"]])
out_summ <- snakemake@output[["summ"]]

nms <- names(X)[grep(".z$", names(X))]
n <- length(nms)

Z <-  X %>%
      select(ends_with(".z")) %>%
      as.matrix()

pvals <- 2*pnorm(-abs(Z))

dot_prds <- expand.grid(n1 = seq(n), n2 = seq(n)) %>%
            filter(n1 <= n2)

prod <- apply(dot_prds, 1, function(x){
                        ix <- which(pvals[,x[1]] > p_thresh &
                                    pvals[,x[2]] > p_thresh)
                        p <- Z[ix,x[1]]*Z[ix,x[2]]
                        s <- sum(p, na.rm=T)
                        nn <- sum(!is.na(p))
                        n1s <- sum(Z[ix,x[1]][!is.na(p)])
                        n2s <- sum(Z[ix,x[2]][!is.na(p)])
                        return(list(s = s, nn=nn, n1s=n1s, n2s=n2s))
                    })
dot_prds$prod <- map(prod, 1) %>% unlist()
dot_prds$n <- map(prod, 2) %>% unlist()
dot_prds$n1s <- map(prod, 3) %>% unlist()
dot_prds$n2s <- map(prod, 4) %>% unlist()
dot_prds$n1 <- nms[dot_prds$n1] %>% unlist()
dot_prds$n2 <- nms[dot_prds$n2] %>% unlist()

saveRDS(dot_prds, file=out_summ)




