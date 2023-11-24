
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
out <- snakemake@out[["out"]]

if(thresh == 1){
  system(paste0("cp ", snakemake@input[["R"]], " ", out))
}else{
  R <- readRDS(snakemake@input[["R"]])
  grps <- corclust(cov2cor(R$R), thresh)
  keep <- sapply(grps, function(x){x[1]})
  newR <- list(R = R$R[keep,keep],
            names = R$names[keep])
  saveRDS(newR, file = out )
}

