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

ss = pd.read_csv(config["input"]["sum_stats"], na_filter=False)

# output options
data_dir = config["out"]["data_dir"] #where the data is
out_dir = config["out"]["output_dir"] #where results will go
formatted_gwas_dir = config["out"]["formatted_gwas_dir"] # where formatted gwas data lives

prefix = config["input"]["label"] + "_"

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

inp = expand(out_dir + prefix + "gfa_{gfas}.ldpruned_{lds}.R_{rs}.final.RDS",
                gfas = gfa_strings,
                lds = ld_strings,
                rs = R_strings)

rule all:
    input: inp

# This produces one data frame per chromosome with columns for snp info
# and columns <study>.z, <study>.ss for z-score and sample size of each snp
rule snp_table_chrom:
    input: files = ss['raw_data_path'], gwas_info = config["input"]["sum_stats"]
    output: out =  temp(data_dir + prefix + "zmat.{chrom}.RDS")
    #params: #nmiss_thresh = config["analysis"]["nmiss_thresh"],
    params: af_thresh = config["analysis"]["af_thresh"],
            sample_size_tol = config["analysis"]["sample_size_tol"]
    wildcard_constraints: chrom = "\d+"
    script: 'R/1_combine_and_format.R'


# LD prune with plink
# LD prune prioritizing snps either by min p-value or  min rank

rule ld_prune_plink:
    input: zmat = data_dir + prefix + "zmat.{chrom}.RDS",
           bfile = config["analysis"]["ldprune"]["ref_path"] + ".bed"
    output: out = temp(data_dir + prefix + "zmat.ldpruned_r2{r2_thresh}_kb{kb}_{p}.{chrom}.RDS")
    params: ref_path = config["analysis"]["ldprune"]["ref_path"],
            pthresh = 1
    wildcard_constraints: chrom = "\d+"
    script: 'R/3_ld_prune_chrom_plink.R'


## Estimate R

# For p-value threshold and ldsc_quick methods, we can compute R
# without ever reading in all of the data.
# For ldsc method, we need to run ldsc for each pair of traits first.

####p-value threshold method

# rule score_summ:
#     input: zmat =  data_dir + prefix + "zmat.ldpruned_r2{r2_thresh}_kb{kb}_{p}.{chrom}.RDS"
#     output: summ = temp(data_dir + prefix + "zmat_summary.ldpruned_r2{r2_thresh}_kb{kb}_{p}.R_pt{pt}.{chrom}.RDS"),
#     wildcard_constraints: chrom = "\d+"
#     script: "R/4_compute_summary.R"
# 
# rule summ_to_cor:
#     input: expand(data_dir + prefix + "zmat_summary.ldpruned_r2{{r2_thresh}}_kb{{kb}}_{{p}}.R_pt{{pt}}.{chrom}.RDS", chrom = range(1, 23))
#     output: out = data_dir + prefix + "R_estimate.ldpruned_r2{r2_thresh}_kb{kb}_{p}.R_pt{pt}.RDS"
#     wildcard_constraints:
#            pt = "[\d.]+"
#     script: "R/4_summary_to_cor.R"

rule R_pt:
  input: Z = expand(data_dir + prefix + "zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23))
  output: out = data_dir + prefix + "R_estimate.ldpruned_r2{r2}_kb{kb}_{p}.R_pt{pt}.RDS"
  wildcard_constraints: pt = "[\d.]+"
  script: "R/4_R_pthresh.R"

###ldsc without updated weights "ldsc_quick"

l2_dir = config["analysis"]["R"]["l2_dir"]
rule score_summ_ldsc:
    input: zmat =  data_dir + prefix + "zmat.{chrom}.RDS",
           l2 = l2_dir + "{chrom}.l2.ldscore.gz"
    output: summ =  temp(data_dir + prefix + "zmat_ldsc_summary.{chrom}.RDS"),
    wildcard_constraints: chrom = "\d+"
    script: "R/4_compute_ldsc_summary.R"

rule summ_to_ldsc_cov:
    input: expand(data_dir + prefix + "zmat_ldsc_summary.{chrom}.RDS", chrom = range(1, 23))
    output: out = data_dir + prefix + "R_estimate.R_ldsc_quick.RDS"
    script: "R/4_ldsc_summ_to_cor.R"

### None
rule none_R:
    input: gwas_info = config["input"]["sum_stats"]
    output: out = data_dir + prefix + "R_estimate.R_none.RDS"
    script: 'R/4_R_none.R'


### Full LDSC compute by pair
# rule ldsc_rg_pair:
#     input: Z = expand(data_dir + prefix + "zmat.{chrom}.RDS", chrom = range(1, 23)),
#            l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23)),
#            m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
#            gwas_info = config["input"]["sum_stats"]
#     output: out =  data_dir + prefix + "ldsc.{name1}___{name2}.RDS"
#     wildcard_constraints:
#       name1 = "(" + "|".join(ss['name']) + ")",
#       name2 = "(" + "|".join(ss['name']) + ")"
#     params: l2_dir = l2_dir
#     script: 'R/4_ldsc_pair.R'
# 
# 
# name_pairs = [(n1, n2) for i1, n1 in enumerate(ss['name']) for i2, n2 in enumerate(ss['name']) if i1 <= i2]
# 
# rule R_ldsc_full:
#     input: data = expand(data_dir + prefix + "ldsc.{np[0]}___{np[1]}.RDS", np = name_pairs),
#            l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23)),
#            gwas_info = config["input"]["sum_stats"]
#     output: out = data_dir + prefix + "R_estimate.R_ldsc.RDS"
#     params: root = data_dir + prefix
#     script: 'R/4_R_ldsc_full.R'
rule R_ldsc_full:
    input: Z = expand(data_dir + prefix + "zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23)), 
           gwas_info = config["input"]["sum_stats"],
           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
    output: out = data_dir + prefix + "R_estimate.R_ldsc.RDS"
    wildcard_constraints: pt = "[\d.]+"
    script: "R/4_R_ldsc_all.R"  

rule cor_clust:
    input: R = data_dir + prefix + "R_estimate.{rstring}.RDS"
    output: out = data_dir + prefix + "R_estimate.{rstring}_cc{cc}.RDS"
    wildcard_constraints:
           cc = "\d+(\.\d+)?"
    script: 'R/4_R_corr_clust.R'
# Run GFA
#

def R_input(wcs):
    global data_dir
    global prefix
    if wcs.Rtype.startswith("pt"):
        return f'{data_dir}{prefix}R_estimate.ldpruned_r2{wcs.r2}_kb{wcs.kb}_{wcs.p}.R_{wcs.Rtype}.RDS'
    else:
        return f'{data_dir}{prefix}R_estimate.R_{wcs.Rtype}.RDS'


rule run_gfa:
    input: Z = expand(data_dir + prefix + "zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23)),
           R = R_input
    output:  out = out_dir + prefix + "gfa_{mode}_gfaseed{fs}_{method}.ldpruned_r2{r2}_kb{kb}_{p}.R_{Rtype}.1.RDS",
    params: params_file = config["analysis"]["gfa_params"],
            max_snps = config["analysis"]["max_snps"]
    wildcard_constraints: fs = "\d+"
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

