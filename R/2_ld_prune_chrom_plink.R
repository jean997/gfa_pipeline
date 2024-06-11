library(dplyr)
library(ieugwasr)

X <- readRDS(snakemake@input[["zmat"]])
r2_thresh <- as.numeric(snakemake@wildcards[["r2_thresh"]])
clump_kb <- snakemake@wildcards[["kb"]]
ref_path  <- snakemake@params[["ref_path"]]
out <- snakemake@output[["out"]]
p <- snakemake@wildcards[["p"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])
is_mvmr <- as.numeric(snakemake@params[["is_mvmr"]])

if(!p %in% c("pvalue", "rank")){
  stop("Unknown prioritization option.\n")
}


Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()

if(is_mvmr == 1){
  Z_hat <- Z_hat[,-1,drop = FALSE]
}

if(p == "pvalue"){
  zmax <- apply(abs(Z_hat), 1, function(x){max(x, na.rm=T)})
  myp <- 2*pnorm(-abs(zmax))

}else if(p == "rank"){
  Z_rank <- apply(Z_hat,2,function(x){rank(x,na.last = "keep")})
  min_rank <- apply(Z_rank, 1, function(x){min(x, na.rm = T)})
  myp <- min_rank/max(min_rank)

}

dat <- data.frame(rsid = X$snp, pval = myp)



dat_clump <- ld_clump(dat = dat,
                     clump_r2 = r2_thresh,
                     clump_p = pthresh,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_path)

ix <- which(X$snp %in% dat_clump$rsid)
X <- X[ix,]

saveRDS(X, file=out)

