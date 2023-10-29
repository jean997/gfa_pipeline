library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(rlang)
library(readr)
library(purrr)
library(stringr)

c <- snakemake@wildcards[["chrom"]]
gwas_info_file <- snakemake@params[["gwas_info"]]
dir <- snakemake@params[["d"]]
nmiss_thresh <- as.numeric(snakemake@params[["nmiss_thresh"]])
af_thesh <- as.numeric(snakemake@params[["af_thresh"]])
out <- snakemake@output[["out"]]


info <- read_csv(gwas_info_file)

fulldat <- map(seq(nrow(info)),   function(i){
                        f <- paste0(dir, info$name[i], ".vcf.bgz")
                        n <- info$name[i]
                        ss <- info$pub_sample_size[i]
                        if(is.na(ss)) ss <- info$pub_cases[i] + info$pub_controls[i]
                        v <- query_chrompos_file(paste0(c, ":1-536870911"), f)
                        pos_name <- as_name(paste0(n, ".pos"))
                        z_name <- as_name(paste0(n, ".z"))
                        ss_name <- as_name(paste0(n, ".ss"))
                        dat <- vcf_to_tibble(v) %>%
                                dplyr::filter(AF > af_thresh & AF < (1-af_thresh)) %>% 
                                dplyr::rename(snp=rsid) %>%
                                dplyr::mutate(Z  = ES/SE) %>%
                                dplyr::select(chr:=seqnames,  snp, REF, ALT, #start, Z, SS)
                                              !!pos_name := start,
                                              !!z_name := Z,
                                              !!ss_name := SS)
                        if(all(is.na(dat[[ss_name]]))) dat[[ss_name]] <- ss
                        return(dat)
                 }) %>%
       purrr::reduce(full_join, by = c("chr", "snp", "REF", "ALT"))

dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
if(length(dup_snps) > 0){
    fulldat <- filter(fulldat, !snp %in% dup_snps)
}



# Save table of how traits are missing each SNP for LD clumping
miss <- fulldat %>%
        select(ends_with(".z")) %>%
        is.na(.) %>%
        rowSums(.)

#nmiss <- data.frame(snp = fulldat$snp, miss = miss)

ix <- which(miss <= nmiss_thresh)

saveRDS(fulldat[ix,], file=out)

#saveRDS(nmiss[ix,], file=nm_out)

