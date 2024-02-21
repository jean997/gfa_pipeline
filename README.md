## Installation Instructions 

These instructions will:
1. Install  Snakemake 
3. Run the GFA pipeline using Snakemake

### Step 1: Installing Snakemake

Follow [instructions here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
to install Snakemake



### Step 2: Install Some necessary R packages

You will need to following R packages. 

- `GFA`: Install using `devtools::install_github("jean997/GFA")`
- `gwasvcf`: Install with `devtools::install_github("mrcieu/gwasvcf")`
- `ieugwasr`: Install with `devtools::install_github("mrcieu/ieugwasr")`
- `stringr`: Install with `install.packages`
- `magrittr`: Install with `install.packages`
- `LaplacesDemon`: Install with `install.packages`. 
- `VariantAnnotation`: Installed from bioconductor using instructions [here](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html).


## Step 3: Using the pipeline

To use the pipeline, create the directory you would like to do your analysis in. Clone the pipeline into that directory using 

```
git clone https://github.com/jean997/gfa_pipeline.git
```

This will create a sub-directory `gfa_pipeline`. We want to move all the files up one level into the working directory. 

```
mv gfa_pipeline/* .
rm -r gfa_pipeline
```

### Create a csv file with all the studies you want to analyze. 

The goal of this section is to create a csv containing all of the studies you wish to analyze with required information.
Your directory should contain a file called `example.csv` which is an example of the file format.
The GWAS data itself can be in one of three formats: 

- A vcf file downloaded from the IEU Open GWAS database
- A flat file (eg csv, tsv etc) with columns for snp, effect and non-effect alleles, effect estimate, and standard error


The pipeline is picky about 
column names but not about the order of columns. The following columns are required. (\*) means that NA values are allowed for any row. (\**) means that NA values 
are only allowed if the data is a vcf file from the IEU open GWAS project or the data are in the Neale format. Since this csv will be read in by R, NA is the preferred way to indicate missing data (as opposed to leaving a blank cell). 

- name: A unique name for the study (you can use `Unique_ID` from the gwas_reference.csv file).
- raw_data_path: Path to the raw data file (`file` in `gwas_reference.csv`).
- neale_format: This column should contain `FALSE` unless the data are GWAS round 2 results from [here](http://www.nealelab.is/uk-biobank). 
- effect_is_or: This value should be "yes" or "no". If yes, this indicates that the value in the beta column is an odds ratio rather than a log odds ratio. 
- pub_sample_size (*): The sample size as stated in the publication or report about the study. 
- snp (**): Column name for rsid
- A1 (**): Column name for effect allele
- A2 (**): Column name for alternate allele
- beta_hat (**): Column name for coefficient estimate
- se (**): Column name for standard error of beta_hat
- p_value (*): Column name for p-value 
- pos (**): Column name for genome position 
- chrom (*): Column name for chromosome 
- sample_size (*): Column name for per-SNP sample size 



Note that the file paths in `raw_data_path` should be absolute paths. 

#### Data from IEU Open GWAS project

These data are the easiest to record because they have a standard format. Include information in the `name` and `raw_data_path` columns. raw_data_path should give the location of the .vcf.gz file. Record information in the `pub_sample_size` column if available. It is a good idea but not required to download the .tbi file and put it in the same directory.

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

This will print the first ten lines of the file. Note that it is important to get the effect and alternate allele columns correct. Sometimes this is easy to figure out because they are called something like "effect_allele" and "other_allele". Sometimes it is more ambiguous and you may need to check any documentation with the file carefully. If columns are called "A1" and "A2" it is almost always the case that A1 is the effect allele and A2 is the non-effect allele. However, it is still good to check this. 

Not all files are usable as they are distributed. For the pipeline to work, we need columns corresponding to rs ID (snp name beginning with "rs"), effect and alternate allele, regression coefficient estimate (beta), and standard error. The vast majority of studies include this information. However, if any of these are missing, you will need to process the file in some way to get it into a usable format. You may not be able to make all formats usable. 

### Editing the config file

At a minimum, you will need to edit the config file so that your csv file is shown on the line for `sum_stats:` in the input section. You may want to modify other options. The config file contains comments describing each option. 


### Running the pipeline

We recommend running the pipeline on a compute cluster. If you are doing this, you may need to edit the `run-snakemake.sh` script to match your cluster architecture. If you are running on a cluster, we like to run 
Snakemake from an interactive session on a compute node. 

To run the pipeline use 

```
./run-snakemake.sh 
```
or to run in the background


```
nohup ./run-snakemake.sh &
```

Using `nohup` and `&` will make it so that the command line returns and all of the Snakemake output is captured in a file called `nohup.out` which is useful if you want to check on the progress. 


## Results

Pipeline results will end up in a folder in your directory called `results` (unless you have changed the `output_dir`) option in `config.yaml`. 
