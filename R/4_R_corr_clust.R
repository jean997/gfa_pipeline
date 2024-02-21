
corclust <- function(R, thresh){

  groups <- list()

  ungrouped <- 1:nrow(R)

  while(length(ungrouped) > 0){
    done <- FALSE
    i <- length(groups) + 1
    groups[[i]] <- c(ungrouped[1])
    ungrouped <- ungrouped[-1]
    grp_len <- length(groups[[i]])
    while(!done){
      Rslice <- abs(R[groups[[i]],,drop = FALSE])
      ix <- which(apply(Rslice, 2, max) > thresh)
      groups[[i]] <- unique(c(groups[[i]], ix))
      if(length(groups[[i]]) == grp_len){
        done <- TRUE
      }else{
        grp_len <- length(groups[[i]])
      }
    }
    ungrouped <- ungrouped[!ungrouped %in% groups[[i]]]
  }
  return(groups)
}





thresh <- as.numeric(snakemake@wildcards[["cc"]])
out <- snakemake@output[["out"]]
cond_num <- as.numeric(snakemake@params[["cond_num"]])-1

  R <- readRDS(snakemake@input[["R"]])
  grps <- corclust(cov2cor(R$R), thresh)
  keep <- sapply(grps, function(x){x[1]})
  newR <- list(R = R$R[keep,keep],
            names = R$names[keep])
  if(!is.null(R$Rg)){newR$Rg <- R$Rg[keep,keep]}
  v <- eigen(newR$R, only.values = TRUE)$values
  if(max(v)/min(v) > cond_num | any(v < 0)){
      newR$R <- GFA::condition(newR$R, cond_num, corr = TRUE)
  }else{
      newR$R <- cov2cor(newR$R)
  }
  newR$groups <- lapply(grps, function(g){R$names[g]})
  saveRDS(newR, file = out)

