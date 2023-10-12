library(stringr)
library(readr)
library(sumstatFactors)
library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(magrittr)

data_file = snakemake@input[["raw_data"]]
snp_name <- snakemake@params[["snp"]]
A1_name <- snakemake@params[["A1"]]
A2_name <- snakemake@params[["A2"]]
beta_hat_name <- snakemake@params[["beta_hat"]]
se_name <- snakemake@params[["se"]]
chrom_name <- snakemake@params[["chrom"]]
pos_name <- snakemake@params[["pos"]]
p_value_name <- snakemake@params[["p_value"]]
effect_is_or <- snakemake@params[["is_or"]]
sample_size_name <- snakemake@params[["sample_size"]]
output_file = snakemake@output[["out"]]


if(str_ends(data_file, "vcf.gz") | str_ends(data_file, "vcf.bgz")){
    #Harmonize strand in vcf to match cause format (A1 = A)
    # Note REF = A2, ALT = A1
    dat <- readVcf(data_file) %>%
           vcf_to_tibble()

    dat <- dat %>%
           rename(A1 = ALT, A2 = REF)
    #remove non snp
    l2 <- str_length(dat$A2)
    l1 <- str_length(dat$A1)
    dat <- dat[l1==1 & l2==1,]

    dat <- sumstatFactors:::remove_ambiguous(dat)
    dat1 <- sumstatFactors:::align_beta(dat, "ES")

    dat <- dat1 %>%
           mutate(AF = case_when(ES == -1*dat$ES ~ 1-AF,
                                 TRUE ~ AF)) %>%
           rename(ALT = A1, REF = A2) %>%
           filter(!is.na(ID))
    out <- dat %$% create_vcf( chrom=seqnames, pos=start,  nea=REF,
                      ea=ALT, snp=ID, ea_af=AF,
                      effect=ES,  se=SE, pval=10^-LP, n=SS, name="a")
    output_file <- str_replace(output_file, ".bgz$", "")
    writeVcf(out, file=output_file, index=TRUE)
}else{
    if(p_value_name != "NA"){
        pstring <- paste0(", `", p_value_name, "`='d'")
    }else{
        pstring <- ""
        p_value_name <- NA
    }
    if(sample_size_name!="NA"){
        sstring <- paste0(", `", sample_size_name, "`='d'")
    }else{
        sstring <- ""
        sample_size_name <- NA
    }
    if(pos_name!="NA"){
        posstring <- paste0(", `", pos_name, "`='d'")
    }else{
        posstring <- ""
        pos_name <- NA
    }

    col_string <- paste0("cols_only(`", snp_name, "`='c', `",
                     A1_name , "`='c', `", A2_name, "`='c', `",
                     beta_hat_name , "`='d', `", se_name, "`='d', `",
                     chrom_name, "`='c' ", posstring,
                     pstring,  sstring, ")")

    X <- read_table2(data_file, col_types = eval(parse(text = col_string)))
    if(effect_is_or == "Yes"){
        X$beta <- log(X[[beta_hat_name]])
        beta_hat <- "beta"
    }

    dat <- gwas_format(X, snp_name, beta_hat_name, se_name, A1_name,
                       A2_name, chrom_name, pos_name,
                       p_value = p_value_name,
                       sample_size = sample_size_name,
                       compute_pval = TRUE)
    #remove non-snp
    l2 <- str_length(dat$A2)
    l2[is.na(l2)] <- 0
    l1 <- str_length(dat$A1)
    l1[is.na(l1)] <- 0
    dat <- dat[l1==1 & l2==1,]

    out <- dat %$% create_vcf( chrom=chrom, pos=pos,  nea=A2,
                      ea=A1, snp=snp,
                      effect=beta_hat,  se=se, pval=p_value, n=sample_size, name="a")
    output_file <- str_replace(output_file, ".bgz$", "")
    writeVcf(out, file=output_file, index=TRUE)
}
