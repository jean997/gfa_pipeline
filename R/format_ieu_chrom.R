require(VariantAnnotation)
require(gwasvcf)
require(dplyr)
require(stringr)
require(readr)
require(GFA)

format_ieu_chrom <- function(file, chrom, af_thresh){
    stopifnot(af_thresh < 0.5)

  dat <- query_chrompos_file(paste0(chrom, ":1-536870911"), file) %>%
    vcf_to_tibble() %>%
    mutate(p_value = 10^{-1*LP})
  if(!all(is.na(dat$AF))){
    dat <- dat %>%
      filter(AF > af_thresh & AF < (1-af_thresh))
  }

  dat <- gwas_format(dat, "ID", "ES", "SE", "ALT",
                       "REF", "seqnames", "start",
                       p_value = "p_value",
                       sample_size = "SS",
                       allele_freq = "AF",
                       compute_pval = TRUE)


  return(dat)
}


format_flat_chrom <- function(file, chrom, af_thresh,
                              snp_name,
                              pos_name,
                              chrom_name,
                              A1_name, A2_name,
                              beta_hat_name,
                              se_name,
                              p_value_name,
                              af_name,
                              sample_size_name,
                              effect_is_or
                              ){

    if(!p_value_name %in% c("", "NA", NA)){
        pstring <- paste0(", `", p_value_name, "`='d'")
    }else{
        pstring <- ""
        p_value_name <- NA
    }
    if(!sample_size_name %in% c("", "NA", NA)){
        sstring <- paste0(", `", sample_size_name, "`='d'")
    }else{
        sstring <- ""
        sample_size_name <- NA
    }
    if(!pos_name %in% c("", "NA", NA)){
        posstring <- paste0(", `", pos_name, "`='d'")
    }else{
        posstring <- ""
        pos_name <- NA
    }
    if(!af_name %in% c("", "NA", NA)){
        afstring <- paste0(", `", af_name, "`='d'")
    }else{
        afstring <- ""
        af_name <- NA
    }


    col_string <- paste0("cols_only(`", snp_name, "`='c', `",
                     A1_name , "`='c', `", A2_name, "`='c', `",
                     beta_hat_name , "`='d', `", se_name, "`='d', `",
                     chrom_name, "`='c' ", posstring,
                     pstring,  sstring, afstring, ")")


    if(str_ends(file, "gz") ){
        h <- read_table(pipe(paste0("gzip -cd ", file, " | head -2")))
        n <- which(names(h) == chrom_name)
        awk_cmd <- paste0("gzip -cd ", file, " | awk '{if ($", n, " == \"", chrom, "\") print $0}' - ")
    }else{
        h <- read_table(pipe(paste0("head -2 ", file)))
        n <- which(names(h) == chrom_name)
        awk_cmd <- paste0("awk '{if ($", n, " == \"", chrom, "\") print $0}' ", file)
    }
    X <- read_table(pipe(awk_cmd), col_types = eval(parse(text = col_string)), col_names = names(h))

    if(!is.na(af_name)){
        if(!all(is.na(X[[af_name]]))){
          ix <- which(X[[af_name]] > af_thresh & X[[af_name]] < (1-af_thresh))
          X <- X[ix,]
        }
    }

    if(effect_is_or){
        X$beta <- log(X[[beta_hat_name]])
        beta_hat <- "beta"
    }

    dat <- gwas_format(X, snp_name, beta_hat_name, se_name,
                       A1 = A1_name,
                       A2 = A2_name,
                       chrom = chrom_name,
                       pos = pos_name,
                       p_value = p_value_name,
                       sample_size = sample_size_name,
                       allele_freq = af_name,
                       compute_pval = TRUE)

    return(dat)
}

