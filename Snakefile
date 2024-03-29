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

localrules: all, check_success, fail, status, final_file, cor_clust
###### Load configuration file
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

# output options
data_dir = config["out"]["data_dir"] #where the data is
out_dir = config["out"]["output_dir"] #where results will go
l2_dir = config["analysis"]["R"]["l2_dir"]


prefix_dict = dict(zip(config["input"]["label"], config["input"]["sum_stats"]))



#ss = pd.read_csv(config["input"]["sum_stats"], na_filter=False)
#prefix = config["input"]["label"] + "_"

gfa_strings = expand("{mode}_gfaseed{s}_{method}",
                     mode = config["analysis"]["mode"],
                     s = config["analysis"]["gfa_seed"],
                     method = config["analysis"]["method"])

ld_strings = expand("r2{r2}_kb{kb}_{p}",
                    r2 = config["analysis"]["ldprune"]["r2_thresh"],
                    kb = config["analysis"]["ldprune"]["clump_kb"],
                    p = config["analysis"]["ldprune"]["ld_prioritization"])

if "pt" in config["analysis"]["R"]["type"]:
    R_type = expand("pt{pt}",
                    pt = config["analysis"]["R"]["pthresh"])
else:
    R_type = []

if "ldsc" in config["analysis"]["R"]["type"]:
    R_type.append("ldsc")

if "ldsc_quick" in config["analysis"]["R"]["type"]:
    R_type.append("ldsc_quick")

if "none" in config["analysis"]["R"]["type"]:
    R_type.append("none")

R_strings = expand("{rt}_cc{cc}",
                   rt = R_type,
                   cc = config["analysis"]["R"]["cor_clust"])

#inp = #expand(out_dir + "{prefix}_gfa_{gfas}.ldpruned_{lds}.R_{rs}.final.RDS",
inp = expand(out_dir + "{prefix}_gls_loadings.{gfas}.ldpruned_{lds}.R_{rs}.Rgcor.RDS",
                prefix = config["input"]["label"],
                gfas = gfa_strings,
                lds = ld_strings,
                rs = R_strings)

rule all:
    input: inp

# This produces one data frame per chromosome with columns for snp info
# and columns <study>.z, <study>.ss for z-score and sample size of each snp
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
    output: out =  data_dir + "{prefix}_zmat.{chrom}.RDS"
    params: af_thresh = config["analysis"]["af_thresh"],
            sample_size_tol = config["analysis"]["sample_size_tol"]
    wildcard_constraints: chrom = "\d+"
    script: 'R/1_combine_and_format.R'


# LD prune with plink
# LD prune prioritizing snps either by min p-value or  min rank

rule ld_prune_plink:
    input: zmat = data_dir + "{prefix}_zmat.{chrom}.RDS",
           bfile = config["analysis"]["ldprune"]["ref_path"] + ".bed"
    output: out = data_dir + "{prefix}_zmat.ldpruned_r2{r2_thresh}_kb{kb}_{p}.{chrom}.RDS"
    params: ref_path = config["analysis"]["ldprune"]["ref_path"],
            pthresh = 1
    wildcard_constraints: chrom = "\d+"
    script: 'R/3_ld_prune_chrom_plink.R'


## Estimate R

# For p-value threshold and ldsc_quick methods, we can compute R
# without ever reading in all of the data.
# For ldsc method, we need to run ldsc for each pair of traits first.

####p-value threshold method

rule R_pt:
  input: Z = expand(data_dir + "{{prefix}}_zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23))
  output: out = data_dir + "{prefix}_R_estimate.ldpruned_r2{r2}_kb{kb}_{p}.R_pt{pt}.RDS"
  wildcard_constraints: pt = "[\d.]+"
  script: "R/4_R_pthresh.R"

###ldsc without updated weights "ldsc_quick"


### None
rule none_R:
    input: gwas_info = info_input
    output: out = data_dir + "{prefix}_R_estimate.R_none.RDS"
    script: 'R/4_R_none.R'


rule R_ldsc_full:
    input: Z = expand(data_dir + "{{prefix}}_zmat.{chrom}.RDS", chrom = range(1, 23)),
           gwas_info = info_input,
           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
    wildcard_constraints: pt = "[\d.]+"
    script: "R/4_R_ldsc_all.R"

rule cor_clust:
    input: R = data_dir + "{prefix}_R_estimate.{rstring}.RDS"
    output: out = data_dir + "{prefix}_R_estimate.{rstring}_cc{cc}.RDS"
    params: cond_num = config["analysis"]["cond_num"]
    wildcard_constraints:
           cc = "\d+(\.\d+)?"
    script: 'R/4_R_corr_clust.R'
# Run GFA
#

def R_input(wcs):
    global data_dir
    if wcs.Rtype.startswith("pt"):
        return f'{data_dir}{wcs.prefix}_R_estimate.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.RDS'
    else:
        return f'{data_dir}{wcs.prefix}_R_estimate.R_{wcs.Rtype}.RDS'


rule run_gfa:
    input: Z = expand(data_dir + "{{prefix}}_zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23)),
           R = R_input
    output:  out = out_dir + "{prefix}_gfa_{mode}_gfaseed{fs}_{method}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.1.RDS",
    params: params_file = config["analysis"]["gfa_params"],
            max_snps = config["analysis"]["max_snps"]
    wildcard_constraints: fs = "\d+"
    script: 'R/5_run_gfa.R'

# This step refits if convergence was not reached
def refit_input(wcs):
    n = int(wcs.n)
    oldn = str(n-1)
    return f'{out_dir}{wcs.prefix}_gfa_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.{oldn}.RDS'

rule refit_gfa:
    input:  inp = refit_input
    output: out = out_dir + "{prefix}_gfa_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.{n}.RDS",
    params: params_file = config["analysis"]["gfa_params"]
    script: 'R/5_refit_gfa.R'

# Below is snakemake machinery for checking convergence and re-running if necessary

max_gfa_tries =  int(config["analysis"]["maxrep"])

def next_input(wcs):
    global max_gfa_tries
    global out_dir

    check_prefix = f'{wcs.prefix}_check_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.'
    check_files = [y for y in os.listdir(out_dir) if y.startswith(check_prefix) ]
    n_tries = len(check_files) + 1

    success_file = f'{out_dir}{wcs.prefix}_success_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.txt'
    fail_file = f'{out_dir}{wcs.prefix}_fail_{wcs.mode}_gfaseed{wcs.fs}_{wcs.pv}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.txt'
    #return fail_file
    if os.path.exists(success_file):
        return success_file
    elif n_tries > max_gfa_tries:
        return fail_file
    else:
        checkpoints.check_success.get(n=n_tries, **wcs)

def final_rds(wcs):
    global out_dir

    check_prefix = f'{wcs.prefix}_check_{wcs.analysis}'
    check_files = [y for y in os.listdir(out_dir) if y.startswith(check_prefix) ]
    n_tries = len(check_files)
    final_gfa_file = f'{out_dir}{wcs.prefix}_gfa_{wcs.analysis}.{n_tries}.RDS'
    return final_gfa_file

checkpoint check_success:
    input: out_dir + "{prefix}_gfa_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.{n}.RDS"
    output: out_check = out_dir + "{prefix}_check_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.{n}.txt"
    params: success_file = out_dir + '{prefix}_success_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.txt'
    wildcard_constraints: n = "\d+"
    script: "R/5_check_gfa.R"

rule fail:
    output: out_dir + '{prefix}_fail_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.txt'
    wildcard_constraints:  n = "\d+"
    params: max_tries = max_gfa_tries
    shell: "echo  Model not converged after {params.max_tries} rounds of {maxiter} iterations. > {output} "

rule status:
    input:  next_input
    output: out = out_dir + "{prefix}_status_{mode}_gfaseed{fs}_{pv}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.txt",
    shell: "cp {input} {output.out}"

rule final_file:
    input: status_file = out_dir + "{prefix}_status_{analysis}.txt",
    output: out_rds = out_dir + "{prefix}_gfa_{analysis}.final.RDS"
    params: rds = final_rds
    shell: "mv {params.rds} {output.out_rds}"

rule gls_loadings_chrom:
    input: z_file = data_dir + "{prefix}_zmat.{chrom}.RDS",
           gfa_file = out_dir + "{prefix}_gfa_{mode}_gfaseed{fs}_{method}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.final.RDS",
           R = R_input
    output: out = out_dir + "{prefix}_gls_loadings.{mode}_gfaseed{fs}_{method}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.{chrom}.RDS"
    script: "R/6_estL_gls.R"


rule R_ldsc_gls:
    input: Z = expand(out_dir + "{{prefix}}_gls_loadings.{{analysis}}.{chrom}.RDS", chrom = range(1, 23)),
           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
    output: out = out_dir + "{prefix}_gls_loadings.{analysis}.Rgcor.RDS"
    script: "R/6_estL_gencor.R"
