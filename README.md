## Installation Instructions for Great Lakes

These instructions will:
1. Install conda
2. Install Snakemake via conda


### Step 1: Installing Conda

Log into Great Lakes. You will now be in your home directory `/home/<username>`. 

From your home directory type:

```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```

This will download a file called `Mambaforge-Linux-x86_64.sh`. You can see the file if you type `ls`.


Type 

```
bash Mambaforge-Linux-x86_64.sh
```

This will install conda. The installer will ask some questions and make you scroll through a use agreement. 
Accept the default installation location. At the end you will get a prompt that asks if you want the installer to initialize Mambaforge3. Say yes to this.

To take effect you need to exit out of Great Lakes and reconnect. Do that now. 

You should now see `(base)` before your terminal prompt. This means that you are now in the "base" conda environment. Later you may 
have other environments but we will just work with base for now. If you ever want to return to the experience of Great Lakes that 
you had before installing conda, type `conda deactivate`. After this command, you will see that `(base)` has disappeared.  To reactivate the base environment, type `conda activate base`. 

Note: This conda installer includes python 3.9.7. This means that python 3.9.7 is installed through the base environment. When the base environment is loaded, typing `python` will now use the python installation that is through conda rather than the python version that is through great lakes. If you want to use the great lakes version, you can deactivate conda but there is probably no reason to do this. Check your python version:
```
python --version
```

(should be 3.9.7)

```
which python
```

Should be `~/miniforge3/bin/python`. 

You now no longer need to use `module load python` to have access to python 3.9.

You will need to reinstall your python modules. We will do that as part of these instructions. 

### Step 2: Install Snakemake

After logging back into Great Lakes and with the base environment active type 

```
mamba install -c conda-forge -c bioconda snakemake
```

You will need to confirm changes by hitting enter at the prompt (or typing Y). This will install snakemake into your base environment. 
It will take a few minutes. Confim that snakmake is now installed by typing 

```
snakemake -h 
```
You should get a lot of help output. If you see "command not found" the something went wrong. 

Last note: In order for the GWAS file tracker to still work, you need to install the wget and pyyaml python modules in your base environment. Type

```
mamba install wget pyyaml
```

That's it. Installation done. 


## Using the pipeline

To use the pipeline, create the directory you would like to do your analysis in. Clone the pipeline into that directory using 

```
git clone git@github.com:jean997/gfa_pipeline.git
```

This will create a sub-directory `gfa_pipeline`. We want to move all the files up one level into the working directory. 

```
mv gfa_pipeline/* .
rm -r gfa_pipeline
```

### Create a csv file with all the studies you want to analyze. 

You should be able to base your csv file on the `gwas_reference.csv` file that tracks all of our data. The pipeline is picky about 
column headers (titles) but not about the order of columns. The following columns are required. Column names below are followed by indicators: (r) means the column must not be NA for any row. (\*) means that NA values are allowed for any row. (\**) means that NA values 
are only allowed if the data is a vcf file from the IEU open GWAS project. Since these data are read in by R, NA is the preferred way to indicate missing data (as opposed to leaving a blank cell). 

- name (r): A unique name for the study (you can use `Unique_ID` from the gwas_reference.csv file).
- raw_data_path (r): Path to the raw data file (`url` in `gwas_reference.csv`).
- pub_sample_size (*): The sample size as stated in the publication or report about the study. 
- snp (**): Column name for rsid
- A1 (**): Column name for effect allele
- A2 (**): Column name for alternate allele
- beta_hat (**): Column name for coefficient estimate
- se (**): Column name for standard error of beta_hat
- p_value (*): Column name for p-value (note optional)
- pos (*): Column name for genome position (note optional)
- chrom (*): Column name for chromosome (note optional)
- sample_size (*): Column name for per-SNP sample size (note optional)
- neale_format (r): This column should contain `FALSE` unless the data are GWAS round 2 results from [here](http://www.nealelab.is/uk-biobank). These files have a different enough format that I had to write a separate parser for them. 
- effect_is_or (r): This value should be "yes" or "no". If yes, this indicates that the value in the beta column is an odds ratio rather than a log odds ratio. 

#### Data from IEU Open GWAS project

These data are the easiest to record because they have a standard format. Include information in the `name` and `raw_data_path` columns. raw_data_path should give the location of the .vcf.gz file. Record information in the `pub_sample_size` column if available. The .tbi file should be present in the same directory. 

For these files the `effect_is_or` column should always be "no" and the `neale_format` column should always be "FALSE". 

All other columns should be NA

#### Data from GWAS round 2 from the Neale Lab

Data downloaded from [here](http://www.nealelab.is/uk-biobank) also have a special format. These can be reported similarly to the IEU files:

Include information in the `name` and `raw_data_path` columns. raw_data_path should give the location of the .tsv.bgz. Record information in the `pub_sample_size` column if available.

For these files the `effect_is_or` column should always be "no" and the `neale_format` column should always be "TRUE". 

All other columns should be NA

#### Data from other sources

Data from any other source need to be in a format that can be read by `read_table2` from the `readr` package. This is essentially any flat file format such as csv, tsv or txt files. It is ok if the file is gzipped or bzipped (ends with `.gz`). 

For these files you need to provide all of the information in columns marked wiht (\**). To find the column names, you may be able to look in a readme that is distributed with the file. The most reliable way to get the information is just to look at the first few lines of the file using `head`. If the file is not gzipped use: 

```
head myfile.tsv
```

If it is gzipped

```
gzip -cd myfile.tsv.gz | head
```

This will print the first ten lines of the file. Note that it is important to get the efffect and alternate allele columns correct. Sometimes this is easy to figure out because they are called something like "effect_allele" and "other_allele". Sometimes it is more ambiguous and you may need to check any documentation with the file carefully. If columns are called "A1" and "A2" it is almost always the case that A1 is the effect allel and A2 is the non-effect allele. However, it is still good to check this. 

Not all files are usable as they are distributed. For the pipeline to work, we need columns corresponding to rs ID (snp name beginning with "rs"), effect and alternate allele, regression coefficient estimate (beta), and standard error. The vast majority of studies include this information. However, if any of these are missing, you will need to process the file in some way to get it into a usable format. You may not be able to make all formats usable. 