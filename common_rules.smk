
## LD pruning options
l2_dir = config["analysis"]["R"]["l2_dir"]
ld_strings = expand("r2{r2}_kb{kb}_{p}",
                    r2 = config["analysis"]["ldprune"]["r2_thresh"],
                    kb = config["analysis"]["ldprune"]["clump_kb"],
                    p = config["analysis"]["ldprune"]["ld_prioritization"])
                    
## R options
if "pt" in config["analysis"]["R"]["type"]:
    R_strings = expand("pt{pt}",
                    pt = config["analysis"]["R"]["pthresh"])
else:
    R_strings = []

if "ldsc" in config["analysis"]["R"]["type"]:
    R_strings.append("ldsc")

if "none" in config["analysis"]["R"]["type"]:
    R_strings.append("none")
    
    
    

# This produces one data frame per chromosome with columns for snp info
# and columns <study>.z, <study>.ss for z-score and sample size of each snp
## rule for getting list of raw data files
def raw_data_input(wcs):
    global prefix_dict
    mycsv = prefix_dict[wcs.prefix]
    ss = pd.read_csv(mycsv, na_filter=False)
    return ss['raw_data_path']

def info_input(wcs):
    global prefix_dict
    return prefix_dict[wcs.prefix]
  
rule snp_table_chrom:
    input: files = raw_data_input, gwas_info = info_input
    output: out =  temp(data_dir + "{prefix}_zmat.{chrom}.RDS")
    params: af_thresh = af_min,
            sample_size_tol = sstol_max
    wildcard_constraints: chrom = r"\d+"
    script: 'R/1_combine_and_format.R'


# LD prune with plink
# LD prune prioritizing snps either by min p-value or  min rank

rule ld_prune_plink:
    input: zmat = data_dir + "{prefix}_zmat.{chrom}.RDS",
           bfile = config["analysis"]["ldprune"]["ref_path"] + ".bed"
    output: out = temp(data_dir + "{prefix}_zmat.ldpruned_r2{r2_thresh}_kb{kb}_{p}.{chrom}.RDS")
    params: ref_path = config["analysis"]["ldprune"]["ref_path"],
            pthresh = 1, 
            is_mvmr = is_mvmr
    wildcard_constraints: chrom = r"\d+"
    script: 'R/2_ld_prune_chrom_plink.R'


## Estimate R

# For p-value threshold and ldsc_quick methods, we can compute R
# without ever reading in all of the data.
# For ldsc method, we need to run ldsc for each pair of traits first.

####p-value threshold method

rule pt_R:
  input: Z = expand(data_dir + "{{prefix}}_zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23))
  output: out = data_dir + "{prefix}_R_estimate.ldpruned_r2{r2}_kb{kb}_{p}.R_pt{pt}.RDS"
  params: cond_num = cond_num
  wildcard_constraints: pt = r"[\d.]+"
  script: "R/3_R_pthresh.R"


### None
rule none_R:
    input: gwas_info = info_input
    output: out = data_dir + "{prefix}_R_estimate.R_none.RDS"
    script: 'R/3_R_none.R'


rule R_ldsc_full:
    input: Z = expand(data_dir + "{{prefix}}_zmat.{chrom}.RDS", chrom = range(1, 23)),
           gwas_info = info_input,
           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
    params: cond_num = cond_num
    wildcard_constraints: pt = r"[\d.]+"
    script: "R/3_R_ldsc_all.R"

