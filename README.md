## Introduction
This repository contains a collection of three Snakemake pipelines for running genetic factor analysis (GFA), multivariable MR (MVMR), and network MR (NESMR). The third of these is for a method that is under development. 

## Installation Instructions 

### Step 1: Installing Snakemake

Follow [instructions here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
to install Snakemake



### Step 2: Install Some necessary R packages

You will need to following R packages for all pipelines. 

- `GFA`: Install using `devtools::install_github("jean997/GFA")`
- `gwasvcf`: Install with `devtools::install_github("mrcieu/gwasvcf")`
- `ieugwasr`: Install with `devtools::install_github("mrcieu/ieugwasr")`
- `dplyr`, `rlang`, `readr`, `stringr`, `purrr`: Install with `install.packages`. 
- `VariantAnnotation`: Installed from bioconductor using instructions [here](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html).

For the MVMR pipeline you may need (if you select the method)

- `TwoSampleMR`
- `MRBEE` 
- `GRAPPLE`
- `esmr`: `devtools::install_github("jean997/esmr")`

For the NESMR pipeline, you will need the `nesmr` branch of the esmr package which can be installed with `devtools::install_github("jean997/esmr", ref = "nesmr")`. The MVMR pipeline will work with either the main branch or the nesmr branch. 

## Step 3: Set up your analysis. 

### 3.1 Create a Directory
I recommend that you create a new directory for each pipeline run. If you want to perform GFA and MVMR on the same data, make two different directories, one for the MVMR analysis and one for the GFA analysis. The reason to do this is that it is easy to forget to change directory names and end up using the wrong data. It is possible to perform multiple pipeline runs in the same directory if you are careful with directory names. 

Once you have the directory you want to use, clone this repository into it using 


```
git clone https://github.com/jean997/gfa_pipeline.git
```

This will create a directory called `gfa_pipeline` in your analysis directory. You can either change to this directory and do your work there, or move everything up one level with `mv gfa_pipeline/* .`. 


### 3.2 Create One or More Study Description Files

All three pipelines rely on the same format for describing GWAS summary statistics. Create a comma separated (csv) file with one row for each trait in your analysis. This is the study description file. Each pipeline can execute multiple analyses, so you may have multiple study description files. If you are running the MVMR pipeline, the outcome trait must be first. Otherwise, the order does not matter. 

The pipeline can handle two types of data formats, vcf files downloaded from the MRC-IEU Open GWAS database, and flat files (e.g. .txt, .tsv, .csv). Flat files can be uncompressed or gzipped (e.g. .csv.gz is ok). 

The study description file should have the following columns:

- `name`: A unique string identifier for the study
- `raw_data_path`: The full path to the summary statistics file
- `pub_sample_size`: Sample size
- `effect_is_or`: Indicator of if effects are reported as odds ratios. Either TRUE or FALSE.
- `chrom`: If using a flat file, the name of the column containing the chromosome
- `pos`: If using a flat file, the name of the column containing the genomic position.
- `snp`: If using a flat file, the name of the column containing rsid. 
- `A1`: If using a flat file, the name of the column containing the effect allele.
- `A2`: If using a flat file, the name of the column containing the other allele.
- `beta_hat`: If using a flat file, the name of the column containing the effect estimate.
- `se`: If using a flat file, the name of the column containing the standard error of the effect estimate.
- `p_value`: If using a flat file, the name of the column containing the p-value.
- `sample_size`: If using a flat file, the name of the column containing the sample size.
- `af`: If using a flat file, the name of the column containing effect allele frequency. 

For studies stored in vcf files, all of the column names can be listed as NA. Studies stored in flat files may be missing p-value, sample size, or allele frequency columns. If these are not present in the file, list NA in that column. The ordering of columns in the study description file does not matter. 

You need to create one study description file for each analysis. For example, if I want to run GFA twice, once with traits 1-10 and a second time with traits 1-12, I need two description files even though there are 10 traits in common. I do not need two copies of the raw GWAS data.

### 3.3 Set up the config file

You will use the config file that corresponds to the pipeline you are running. If you are running GFA, edit `config_gfa.yaml`. If you are running MVMR, edit `config_mvmr.yaml`. If you are running NESMR, edit `config_nesmr.yaml`. 

Each config file contains three main sections: `input` describing input files, `analysis` describing analysis options and `otuput` specifying output directories. 

#### 3.3.1 The `input` section
For GFA and MVMR, for each analysis include an entry of the form `analysis_name: "description_file.csv"` in the input section. For example, if you want to run two analyses corresponding to two description files, your input section may look like 

```
input:
    test1: "test_data/test_traits1.csv"
    test2: "test_data/test_traits2.csv"
```
The names `test1` and `test2` are arbitrary. You can use something informative about your analysis. You can have as many entries here as you want. 

For NESMR, you need to provide both the study description file and a list of template files. For NESMR, your input section will look something like 

```
input:
    test1:
      trait_csv: "test_data/test_traits1.csv"
      templates: ["test_data/trait1_template1.RDS", "test_data/trait1_template2.RDS"]
    test2:
      trait_csv: "test_data/test_traits2.csv"
      templates: ["test_data/trait2_template1.RDS", "test_data/trait2_template2.RDS", "test_data/trait2_template3.RDS"]
```
Each template file includes one direct effect template stored as an RDS. The pipeline will run one analysis for each template and study combination. The order of the traits in the template matrix corresponds to the order of traits in description file. 

#### 3.3.2 The `analysis` section

All three pipelines have two common sections, `ldprune` describing how the data are LD-pruned and `R` describing how to estimate the nuisance correlation matrix. They also have unique analysis-specific sections. A description of all options can be found in comments in the config file. If an option is described as "list ok", this means that you can provide a list of different parameters. For example, you could use different LD-pruning thresholds by using `r2: [0.01, 0.001]`. This will cause the pipeline to run one analysis for each combination of listed parameters. 

#### 3.3.3 The `output` section

All three pipelines produce two types of output, "data" and "results". Specify directories to store these in. 

## Step 4: Run Your Pipeline

### Note about temporary files

Each of the pipelines will create a full copy of all of the GWAS data in the formatting step. In order to prevent accumulating many copies of data, the pipeline treats these as temporary files and deletes them when the analysis is completed. However, if you are planning to add on to an analysis (for example for NESMR if you think you will add more templates, or for MVMR you think you might like to add more analysis options), it can be useful to keep the data files so that you don't have to recreate them when you run with new options. To prevent deleting temporary data files add the `--notemp` flag to the snakemake command.

### 4.1 Run the pipeline at the command line
To run the pipeline at the command line, first check that all is well with your config files by running a dry run using 

```
snakemake -n -s Snakefile_{pipeline}
```
Fill in the correct ending for your desired pipeline in place of `{pipeline}` above. This should run without error and report the steps the pipeline will run. If you get an error, there may be something wrong with your configuration file. 

Next run the pipeline

```
snakemake -c1 -s Snakefile_{pipeline}
```

The flag `-c1` specifies to use one core. You can use more cores and run more jobs in parallel using `-c n` for n cores. 

### 4.2 Run the pipeline on a compute cluster

Running at the command line is ok if you are using a small number of traits, but we generally recommend running the pipelines on a compute cluster. 

If you are running on a compute cluster, it is still a good idea to do a dry run at the command line as above to make sure everything is set up to run. 

Next, take a look at the `cluster.yaml` file. This file specifies the cluster resources alotted to each step. The values are set at reasonable numbers that work for many analyses. However, if you find that you are getting errors related to running out of time or memory, you may need to change these. You may be able to use less memory for smaller analyses.  Next, take a look at the file `run-snakemake-{pipeline}.sh`. This file contains a command for executing snakemake so that it will submit jobs to your cluster. If your cluster doesn't use slurm, you will need to edit the argument to `--cluster`. Take a look at the Snakemake documentation for more details. If your cluster does use slurm, this command will probably work for you, but you will probably need to either change or delete the `--account` flag. 

To run the pipeline use, 
```
./run-snakemake-{pipeline}.sh 
```
or to run in the background

```
nohup ./run-snakemake-{pipeline}.sh &
```

Using `nohup` and `&` will make it so that the command line returns and all of the Snakemake output is captured in a file called `nohup.out`.




## Results

Pipeline results will end up in the folder specified in the output options of your config file. 

For the GFA pipeline, you will see a lot of files. For each analysis, the main GFA results are in a file ending in `.final.RDS`. You will also see files for GLS weights calculated by chromosome and the genetic correlation of the estimated factors. You will also see some text files. These are used by the pipeline to check convergence and can be ignored. 

For the MVMR and NESMR pipeline, you will see one output file per analysis. The file names should be fairly self explanatory. For example the file `test1_ivw_5e-08.ldpruned_r20.01_kb1000_pvalue.R_ldsc.RDS` is the analysis of the traits in the `test1` description file using IVW-regression with a p-value threshold of `5e-8`. The rest of the filename describes the LD-pruning and nuisance parameter estimation options.
