# Snakemake pipeline for analyzing gwas summary statistic data using GFA
#
#
# LICENSE: CC0. Do what you want with the code, but it has no guarantees.
#          https://creativecommons.org/share-your-work/public-domain/cc0/
#
#
# source activate cause_large
#
# ./run_snakemake.sh
#
# don't forget to update cluster.yaml


import pandas as pd
import random
import string
from snakemake.utils import validate

localrules: all, summ_to_cor, check_success, fail, status, final_file, scale_pve
###### Load configuration file
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

ss = pd.read_csv(config["input"]["sum_stats"], na_filter=False)

# output options
data_dir = config["out"]["data_dir"] #where the data is
out_dir = config["out"]["output_dir"] #where results will go
formatted_gwas_dir = config["out"]["formatted_gwas_dir"] # where formatted gwas data lives

prefix = config["input"]["label"] + "_"

gfa_strings = expand("{mode}_gfaseed{s}_{pv}",
                     mode = config["analysis"]["mode"],
                     s = config["analysis"]["gfa_seed"], 
                     pv = config["analysis"]["pvalue_threshold"])

ld_strings = expand("r2{r2}_kb{kb}_{p}",
                    r2 = config["analysis"]["ldprune"]["r2_thresh"],
                    kb = config["analysis"]["ldprune"]["clump_kb"],
                    p = config["analysis"]["ldprune"]["ld_prioritization"])

if "pt" in config["analysis"]["R"]["type"]:
    R_strings = expand("pt{pt}", pt = config["analysis"]["R"]["pthresh"])
else:
    R_strings = []

if "ldsc" in config["analysis"]["R"]["type"]:
    R_strings.append("ldsc")

if "ldsc_full" in config["analysis"]["R"]["type"]:
    R_strings.append("ldsc_full")

if "none" in config["analysis"]["R"]["type"]:
    R_strings.append("none")


inp = expand(out_dir + prefix + "gfa_{gfas}.ldpruned_{lds}.R_{rs}.scaled.RDS",
                gfas = gfa_strings,
                lds = ld_strings,
                rs = R_strings)

rule all:
    input: inp
    #input: "data/bc_zmat.22.RDS"

# This produces one data frame per chromosome with columns for snp info
# and columns <study>.z, <study>.ss for z-score and sample size of each snp
# The nmiss file has two columns, one for snp one for number of missing studies.
rule snp_table_chrom:
    input: files = ss['raw_data_path'], gwas_info = config["input"]["sum_stats"]
    output: out =  data_dir + prefix + "zmat.{chrom}.RDS"
    params: nmiss_thresh = config["analysis"]["nmiss_thresh"],
            af_thresh = config["analysis"]["af_thresh"]
    wildcard_constraints: chrom = "\d+"
    script: 'R/1_combine_and_format.R'


# LD prune with plink
# LD prune prioritizing snps with less missingness
#ld_pt = max(config["analysis"]["pvalue_threshold"])
ld_pt = 1

rule ld_prune_plink:
    input: zmat = data_dir + prefix + "zmat.{chrom}.RDS",
           bfile = config["analysis"]["ldprune"]["ref_path"] + ".bed"
    output: out = data_dir + prefix + "zmat.ldpruned_r2{r2_thresh}_kb{kb}_{p}.{chrom}.RDS"
    params: ref_path = config["analysis"]["ldprune"]["ref_path"], 
            pthresh = ld_pt
    wildcard_constraints: chrom = "\d+"
    script: 'R/3_ld_prune_chrom_plink.R'


## Estimate R
# We can compute the p-thresholded R matrix without ever reading in all of the data
# We first compute per chromosome summaries and the compute R
# There is an LD score version of this but it does not update the weights and so it has higher variance
# than the published version
# Finally there is the full LD score version

####p-value threshold method

rule score_summ:
    input: zmat =  data_dir + prefix + "zmat.ldpruned_r2{r2_thresh}_kb{kb}_{p}.{chrom}.RDS"
    output: summ =  data_dir + prefix + "zmat_summary.ldpruned_r2{r2_thresh}_kb{kb}_{p}.R_pt{pt}.{chrom}.RDS",
    wildcard_constraints: chrom = "\d+"
    script: "R/4_compute_summary.R"

rule summ_to_cor:
    input: expand(data_dir + prefix + "zmat_summary.ldpruned_r2{{r2_thresh}}_kb{{kb}}_{{p}}.R_pt{{pt}}.{chrom}.RDS", chrom = range(1, 23))
    output: out = data_dir + prefix + "R_estimate.ldpruned_r2{r2_thresh}_kb{kb}_{p}.R_pt{pt}.RDS"
    script: "R/4_summary_to_cor.R"

###ldsc without updated weights

l2_dir = config["analysis"]["R"]["l2_dir"]
rule score_summ_ldsc:
    input: zmat =  data_dir + prefix + "zmat.{chrom}.RDS",
           l2 = l2_dir + "{chrom}.l2.ldscore.gz"
    output: summ =  data_dir + prefix + "zmat_ldsc_summary.{chrom}.RDS",
    wildcard_constraints: chrom = "\d+"
    script: "R/4_compute_ldsc_summary.R"

rule summ_to_ldsc_cov:
    input: expand(data_dir + prefix + "zmat_ldsc_summary.{chrom}.RDS", chrom = range(1, 23))
    output: out = data_dir + prefix + "R_estimate.R_ldsc.RDS"
    script: "R/4_ldsc_summ_to_cor.R"

### None
rule none_R:
    output: out = data_dir + "none_R.txt"
    shell: "touch {output.out}"


### Full LDSC compute by pair
# M doesn't matter so this could be modified to leave it out.
rule ldsc_rg_pair:
    input: f1 = formatted_gwas_dir + "{name1}.vcf.bgz",
           f2 = formatted_gwas_dir + "{name2}.vcf.bgz",
           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23)),
           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23))
    output: out =  data_dir + prefix + "ldsc.{name1}.{name2}.RDS"
    params: gwas_info = config["input"]["sum_stats"],
            l2_dir = l2_dir
    script: 'R/4_ldsc_pair.R'


name_pairs = [(n1, n2) for i1, n1 in enumerate(ss['name']) for i2, n2 in enumerate(ss['name']) if i1 <= i2]

rule R_ldsc_full:
    input: data = expand(data_dir + prefix + "ldsc.{np[0]}.{np[1]}.RDS", np = name_pairs),
           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
    output: out = data_dir + prefix + "R_estimate.R_ldsc_full.RDS"
    params: gwas_info = config["input"]["sum_stats"], root = data_dir + prefix
    script: 'R/4_create_R_ldsc_full.R'

# Run GFA
#

def R_input(wcs):
    global data_dir
    global prefix
    if wcs.Rtype == "ldsc":
        return f'{data_dir}{prefix}R_estimate.R_ldsc.RDS'
    elif wcs.Rtype == "none":
        return f'{data_dir}none_R.txt'
    elif wcs.Rtype == "ldsc_full":
        return f'{data_dir}{prefix}R_estimate.R_ldsc_full.RDS'
    else:
        return f'{data_dir}{prefix}R_estimate.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.RDS'


rule run_gfa:
    input: Z = expand(data_dir + prefix + "zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23)),
           R = R_input
    output:  out = out_dir + prefix + "gfa_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.1.RDS",
    params: params_file = config["analysis"]["gfa_params"],
            max_snps = config["analysis"]["max_snps"]
    script: 'R/5_run_gfa.R'

# This step refits if convergence was not reached
def refit_input(wcs):
    n = int(wcs.n)
    oldn = str(n-1)
    return f'{out_dir}{prefix}gfa_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.{oldn}.RDS'

rule refit_gfa:
    input:  refit_input
    output: out = out_dir + prefix + "gfa_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.{n}.RDS",
    params: params_file = config["analysis"]["gfa_params"]
    script: 'R/5_refit_gfa.R'

# Below is snakemake machinery for checking convergence and re-running if necessary

max_gfa_tries =  int(config["analysis"]["maxrep"])

def next_input(wcs):
    global max_gfa_tries
    global out_dir
    global prefix

    check_prefix = f'{prefix}check_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.'
    check_files = [y for y in os.listdir(out_dir) if y.startswith(check_prefix) ]
    n_tries = len(check_files) + 1

    success_file = f'{out_dir}{prefix}success_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.txt'
    fail_file = f'{out_dir}{prefix}fail_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.txt'
    #return fail_file
    if os.path.exists(success_file):
        return success_file
    elif n_tries > max_gfa_tries:
        return fail_file
    else:
        checkpoints.check_success.get(n=n_tries, **wcs)

def final_rds(wcs):
    global out_dir
    global prefix

    check_prefix = f'{prefix}check_{wcs.analysis}'
    check_files = [y for y in os.listdir(out_dir) if y.startswith(check_prefix) ]
    n_tries = len(check_files)
    final_gfa_file = f'{out_dir}{prefix}gfa_{wcs.analysis}.{n_tries}.RDS'
    return final_gfa_file

checkpoint check_success:
    input: out_dir + prefix + "gfa_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.{n}.RDS"
    output: out_check = out_dir + prefix + "check_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.{n}.txt"
    params: success_file = out_dir + prefix + 'success_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.txt'
    wildcard_constraints: n = "\d+"
    script: "R/5_check_gfa.R"

rule fail:
    output: out_dir + prefix + 'fail_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.txt'
    wildcard_constraints:  n = "\d+"
    params: max_tries = max_gfa_tries
    shell: "echo  Model not converged after {params.max_tries} rounds of {maxiter} iterations. > {output} "

rule status:
    input:  next_input 
    output: out = out_dir + prefix + "status_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.txt",
    shell: "cp {input} {output.out}"

rule final_file:
    input: status_file = out_dir + prefix + "status_{analysis}.txt",
    output: out_rds = out_dir + prefix + "gfa_{analysis}.final.RDS"
    params: rds = final_rds 
    shell: "mv {params.rds} {output.out_rds}"

rule scale_pve:
    input: in_rds =  out_dir + prefix + "gfa_{mode}_gfaseed{str}.final.RDS" 
    output: out = out_dir + prefix + "gfa_{mode}_gfaseed{str}.scaled.RDS"
    params: gwas_info = config["input"]["sum_stats"]
    script: "R/6_scale_pve.R"

#rule estimate_L:
#    input: inp = out_dir + prefix + "fit_{key}.RDS", zmat = data_dir + prefix + "zmat.{chrom}.RDS"
#    output: out = out_dir + prefix + "estL_{key}.{chrom}.RDS"
#    shell: "Rscript R/7_estimate_L.R {input.inp} {input.zmat} {output.out}"
