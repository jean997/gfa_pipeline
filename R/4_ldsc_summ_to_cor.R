library(dplyr)
library(purrr)
library(LaplacesDemon)

args <- commandArgs(trailingOnly=TRUE)
out <-  args[1]
files <-  args[-1]



df <- map(files, function(f){
      readRDS(f) 
      }) %>%
  reduce(bind_rows)  %>%
      group_by(n1, n2) %>%
      summarize(xtx = sum(xtx), 
                xty = sum(xty), 
                m = sum(m), 
                ysum  = sum(ysum), 
                xsum = sum(xsum)) %>% 
      mutate(b1 = (xty - (xsum*ysum/m))/(xtx - ((xsum^2)/m)), 
             b0 = (ysum/m) - b1*(xsum/m))


# we want a symmetric matrix so we need to replicate some rows of df
df_copy <- filter(df, n1 != n2) %>%
           rename(n1c = n2, n2c = n1) %>%
           rename(n1 = n1c, n2 = n2c)

cov_mat <- bind_rows(df, df_copy)  %>%
           select(n1, n2, b0) %>%
           reshape2::dcast(n1 ~ n2)

nms <- cov_mat$n1
R <- as.matrix(cov_mat[,-1])   %>% as.positive.definite()
eS <- eigen(R)

ret <- list(R = R, names = nms, eS = eS)

saveRDS(ret, file=out)


