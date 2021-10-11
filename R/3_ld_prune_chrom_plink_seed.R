library(tidyverse)
library(ieugwasr)

args <- commandArgs(trailingOnly=TRUE)
nmiss <- readRDS(args[1])
r2_thresh <- as.numeric(args[2])
clump_kb <- args[3]
ref_path  <- args[4]
out_nm <- args[5]
seed <- args[6]

nmiss <- nmiss %>% 
         filter(miss == 0 ) %>%
         rename(rsid = snp)


if(seed == 0){
    nmiss <- nmiss %>% 
         mutate(pval = 0.5)
}else{
    set.seed(seed)
    nmiss$pval <- runif(n = nrow(nmiss))
}

nm_clump <- ld_clump(dat = nmiss, 
                     clump_r2 = r2_thresh, 
                     clump_p = 1, 
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(), 
                     bfile = ref_path)

nm_clump <- rename(nm_clump, snp=rsid)

saveRDS(nm_clump, file=out_nm)

