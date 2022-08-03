# Snakemake pipeline for analyzing gwas summary statistic data using flash
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

localrules: all, summ_to_cor, check_success, fail, final_flash
###### Load configuration file
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

ss = pd.read_csv(config["input"]["sum_stats"], na_filter=False)

# output options
data_dir = config["out"]["data_dir"] #where the data is
out_dir = config["out"]["output_dir"] #where results will go
formatted_gwas_dir = config["out"]["formatted_gwas_dir"]

prefix = config["input"]["label"]

gfa_strings = expand("{mode}_gfaseed{s}",
                     mode = config["analysis"]["mode"],
                     s = config["analysis"]["gfa_seed"])

ld_strings = expand("r2{r2}_kb{kb}_seed{s}",
                    r2 = config["analysis"]["ldprune"]["r2_thresh"],
                    kb = config["analysis"]["ldprune"]["clump_kb"],
                    s = config["analysis"]["ldprune"]["ld_seed"])

if "pt" in config["analysis"]["R"]["type"]:
    R_strings = expand("pt{pt}", pt = config["analysis"]["R"]["pthresh"])
else:
    R_strings = []
if "ldsc" in config["analysis"]["R"]["type"]:
    R_strings = R_strings.append("ldsc")

if "none" in config["analysis"]["R"]["type"]:
    R_strings = R_strings.append("none")


inp = expand(out_dir + prefix + "status_{gfas}.ldpruned_{lds}.R_{rs}.txt",
                gfas = gfa_strings,
                lds = ld_strings,
                rs = R_strings)

#est_L = config["analysis"]["est_L"]
#if est_L:
#    est_L_key = config["analysis"]["est_L_key"]
#    inp = expand(out_dir + prefix + "estL_" + est_L_key + ".{chr}.RDS", chr = range(1, 23))


rule all:
    input: inp

# convert data to standardized format and harmonize
rule format:
    input: raw_data  = lambda wildcards: ss[ss['name'] == wildcards.name]['raw_data_path'].tolist()[0]
    output: out = formatted_gwas_dir + "{name}.vcf.bgz"
    params:
        snp = lambda wildcards: ss[ss['name'] == wildcards.name]['snp'].tolist()[0],
        A1 = lambda wildcards: ss[ss['name'] == wildcards.name]['A1'].tolist()[0],
        A2 = lambda wildcards: ss[ss['name'] == wildcards.name]['A2'].tolist()[0],
        beta_hat = lambda wildcards: ss[ss['name'] == wildcards.name]['beta_hat'].tolist()[0],
        se = lambda wildcards: ss[ss['name'] == wildcards.name]['se'].tolist()[0],
        chrom = lambda wildcards: ss[ss['name'] == wildcards.name]['chrom'].tolist()[0],
        pos = lambda wildcards: ss[ss['name'] == wildcards.name]['pos'].tolist()[0],
        p_value = lambda wildcards: ss[ss['name'] == wildcards.name]['p_value'].tolist()[0],
        sample_size = lambda wildcards: ss[ss['name'] == wildcards.name]['sample_size'].tolist()[0],
        is_or = lambda wildcards: ss[ss['name'] ==wildcards.name]['effect_is_or'].tolist()[0],
        neale_format = lambda wildcards: ss[ss['name'] ==wildcards.name]['neale_format'].tolist()[0],
        neale_var_ref = config["input"]["neale_var_ref"]
    shell: 'Rscript R/1_format_data_withukb.R {input.raw_data} {params.snp} {params.A1} \
            {params.A2} {params.beta_hat} {params.se} {params.chrom} {params.pos} {params.p_value} \
            {params.is_or} {params.sample_size} {output.out} {params.neale_format} {params.neale_var_ref}'

# This produces one data frame per chromosome with columns for snp info
# and columns <study>.z, <study>.ss for z-score and sample size of each snp
# The nmiss file has two columns, one for snp one for number of missing studies.
rule snp_table_chrom:
    input: files = expand(formatted_gwas_dir + "{name}.vcf.bgz", name = ss['name'])
    output: out =  data_dir + prefix + "zmat.{chrom}.RDS",
            nmiss = data_dir + prefix + "nmiss.{chrom}.RDS"
    params: gwas_info = config["input"]["sum_stats"], d = formatted_gwas_dir, nmiss_thresh = config["analysis"]["nmiss_thresh"]
    wildcard_constraints: chrom = "\d+"
    shell: 'Rscript R/2_combine_data_vcf.R {wildcards.chrom} {params.gwas_info} \
           {params.d} {params.nmiss_thresh} {output.out} {output.nmiss}'


# LD prune with plink
# LD prune prioritizing snps with less missingness
rule ld_prune_plink:
    input: nmiss = data_dir + prefix + "nmiss.{chrom}.RDS",
           bfile = config["analysis"]["ldprune"]["ref_path"] + ".bed"
    output: out = data_dir + prefix + "nmiss.ldpruned_r2{r2_thresh}_kb{kb}_seed{s}.{chrom}.RDS"
    params: ref_path = config["analysis"]["ldprune"]["ref_path"]
    wildcard_constraints: chrom = "\d+"
    shell:   'Rscript R/3_ld_prune_chrom_plink_seed.R {input.nmiss}   \
                   {wildcards.r2_thresh} {wildcards.kb} {params.ref_path} {output.out} \
                   {wildcards.s}'

# After this step, we don't need to keep the full data
rule subset_zmat:
    input: normbeta =  data_dir + prefix + "zmat.{chrom}.RDS",
           keep = data_dir + prefix + "nmiss.ldpruned_r2{r2_thresh}_kb{kb}_seed{s}.{chrom}.RDS"
    output: normbeta = data_dir + prefix + "zmat.ldpruned_r2{r2_thresh}_kb{kb}_seed{s}.{chrom}.RDS"
    wildcard_constraints: chrom = "\d+"
    shell: 'Rscript R/3_subset_zmat.R {input.normbeta} {input.keep} {output.normbeta}'

## Estimate R
rule score_summ:
    input: normbeta =  data_dir + prefix + "zmat.ldpruned_r2{r2_thresh}_kb{kb}_seed{s}.{chrom}.RDS"
    output: summ =  data_dir + prefix + "zmat_summary.ldpruned_r2{r2_thresh}_kb{kb}_seed{s}.R_pt{pt}.{chrom}.RDS",
    wildcard_constraints: chrom = "\d+"
    shell: 'Rscript R/3_compute_summary.R {input.normbeta} {wildcards.pt} {output.summ}'

l2_dir = config["analysis"]["R"]["l2_dir"]
rule score_summ_ldsc:
    input: normbeta =  data_dir + prefix + "zmat.{chrom}.RDS",
           l2 = l2_dir + "{chrom}.l2.ldscore.gz"
    output: summ =  data_dir + prefix + "zmat_ldsc_summary.{chrom}.RDS",
    wildcard_constraints: chrom = "\d+"
    shell: 'Rscript R/3_compute_ldsc_summary.R {input.normbeta} {input.l2} {output.summ}'

# Compute R from the per chromosome summaries
rule summ_to_cor:
    input: expand(data_dir + prefix + "zmat_summary.ldpruned_r2{{r2_thresh}}_kb{{kb}}_seed{{s}}.R_pt{{pt}}.{chrom}.RDS", chrom = range(1, 23))
    output: out = data_dir + prefix + "R_estimate.ldpruned_r2{r2_thresh}_kb{kb}_seed{s}.R_pt{pt}.RDS"
    shell: "Rscript R/4_summary_to_cor.R  {output.out} {input}"

rule summ_to_ldsc_cov:
    input: expand(data_dir + prefix + "zmat_ldsc_summary.{chrom}.RDS", chrom = range(1, 23))
    output: out = data_dir + prefix + "R_estimate.R_ldsc.RDS"
    shell: "Rscript R/4_ldsc_summ_to_cor.R  {output.out} {input}"

rule none_R:
    output: out = data_dir + "none_R.txt"
    shell: "touch {output.out}"

# Run GFA
#

def R_input(wcs):
    global data_dir
    global prefix
    if wcs.Rtype == "ldsc":
        return f'{data_dir}{prefix}R_estimate.R_ldsc.RDS'
    else if wcs.Rtype == "none":
        return f'{data_dir}none_R.txt'
    else:
        return f'{data_dir}{prefix}R_estimate.ldpruned_r2{wcs.r2}_kb{wcs.kb}_seed{wcs.ls}.R_{wcs.Rtype}.RDS'


rule run_gfa:
    input: NB = expand(data_dir + prefix + "zmat.ldpruned_r2{{r2}}_kb{{kb}}_seed{{ls}}.{chrom}.RDS", chrom = range(1, 23)),
           R = R_input
    output:  out = out_dir + prefix + "gfa_{mode}_gfaseed{fs}.ldpruned_r2{r2}_kb{kb}_seed{ls}.R_{Rtype}.1.RDS",
    params: params_file = config["analysis"]["gfa_params"]
    shell: 'Rscript R/5_run_gfa.R {output.out} {wildcards.mode} {input.R}  \
            {params.params_file} {wildcards.fs} \
            {input.NB}'

# This step refits if convergence was not reached
def refit_input(wcs):
    n = int(wcs.n)
    oldn = str(n-1)
    return f'{out_dir}{prefix}gfa_{wcs.mode}_gfaseed{wcs.fs}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_seed{wcs.ls}.R_{wcs.Rtype}.{oldn}.RDS'

rule refit_gfa:
    input:  refit_input
    output: out = out_dir + prefix + "gfa_{mode}_gfaseed{fs}.ldpruned_r2{r2}_kb{kb}_seed{ls}.R_{Rtype}.{n}.RDS",
    params: params_file = config["analysis"]["gfa_params"]
    shell: 'Rscript R/6_refit_gfa.R {input} {output.out}  {params.param_file}'

# Below is snakemake machinery for checking convergence and re-running if necessary

max_gfa_tries =  int(config["analysis"]["maxrep"])

def next_input(wcs):
    global max_gfa_tries
    global out_dir
    global prefix

    check_prefix = f'{prefix}check_{wcs.mode}_gfaseed{wcs.fs}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_seed{wcs.ls}.R_{wcs.Rtype}.'
    check_files = [y for y in os.listdir('results') if y.startswith(check_prefix) ]
    flash_tries = len(check_files) + 1

    success_file = f'{out_dir}{prefix}success_{mode}_gfaseed{wcs.fs}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_seed{wcs.ls}.R_{wcs.Rtype}.txt'
    fail_file = f'{out_dir}{prefix}fail_{mode}_gfaseed{wcs.fs}.ldpruned_r2{wcs.r2}_kb{wcs.kb}_seed{wcs.ls}.R_{wcs.Rtype}.txt'
    #return fail_file
    if os.path.exists(success_file):
        return success_file
    elif flash_tries > max_gfa_tries:
        return fail_file
    else:
        checkpoints.check_success.get(n=flash_tries, **wcs)

checkpoint check_success:
    input: out_dir + prefix + "fit_{m}_{k}_{ms}_{fm}_{maxiter}_seed{fs}.ldpruned_r2{r2}_kb{kb}_seed{ls}.R_{Rtype}.{n}.RDS"
    output: out = out_dir + prefix + "check_{m}_{k}_{ms}_{fm}_{maxiter}_seed{fs}.ldpruned_r2{r2}_kb{kb}_seed{ls}.R_{Rtype}.{n}.txt",
    params: success_file = out_dir + prefix + 'success_{m}_{k}_{ms}_{fm}_{maxiter}_seed{fs}.ldpruned_r2{r2}_kb{kb}_seed{ls}.R_{Rtype}.txt'
    wildcard_constraints: n = "\d+"
    shell: "Rscript R/99_check_flash.R {input} {wildcards.n} {output.out} {params.success_file}"

rule fail:
    output: out_dir + prefix + 'fail_{mode}_gfaseed{fs}.ldpruned_r2{r2}_kb{kb}_seed{ls}.R_{Rtype}.txt'
    wildcard_constraints:  n = "\d+"
    params: max_tries = max_gfa_tries
    shell: "echo  Model not converged after {params.max_tries} rounds of {maxiter} iterations. > {output} "

rule final_flash:
    input: next_input
    output: out = out_dir + prefix + "status_{mode}_gfaseed{fs}.ldpruned_r2{r2}_kb{kb}_seed{ls}.R_{Rtype}.txt",
    shell: "cp {input} {output}"

rule estimate_L:
    input: inp = out_dir + prefix + "fit_{key}.RDS", zmat = data_dir + prefix + "zmat.{chrom}.RDS"
    output: out = out_dir + prefix + "estL_{key}.{chrom}.RDS"
    shell: "Rscript R/7_estimate_L.R {input.inp} {input.zmat} {output.out}"
