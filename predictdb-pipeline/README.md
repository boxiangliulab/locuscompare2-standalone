# PredictDB-Pipline

Training prediction models and putting them in a format that can be used to run the predixcan.
Original pipline is [here](https://github.com/hakyimlab/PredictDB-Tutorial)

## 1. Setup Environment

1) Required R packages: tidyverse, dplyr, RSQLite, glmnet, reshape2
2) Peer factors covariates generation tool: [peer](https://github.com/hakyimlab/peer)
3) Python packages: pandas

## 2. Input file

+ QTL genotype in [vcf](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format

  ID column in vcf must NOT be missing if you want to identify a varint by rsid in result db and covariances with parameter --use_rsid_in_covariances

+ phenotype file in [bed](https://qtltools.github.io/qtltools/) format.

  The format is a custom UCSC [bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format;
  bed format should have 6+N columns, 1 of the 6 columns must be 'pid' (which means gene_id), N columns are expression values for genes.
  if your bed format is different, use --exp_gene_id_col_name to specify gene_id column name and --exp_rest_non_individual_col_names to specify rest column names(except gene_id and expression value columns)

+ genecode file

  this file can be downloaded from [here](https://www.gencodegenes.org/human/releases.html), pick the version that matching your expression file.


## 3. Run

+ Run the pipline from command line, this is time consuming, it can take hours up to days to complete depending on you data size.

```shell
python3 predictdb-pipeline/pipeline.py --genotype genotpe_vcf_path --expression gene_exp_bed_path --genecode genecode_path --output_dir output_directory
```

## 4. Output

Results model db is placed in [output_directory]/db, covariances is in [output_directory]/covariances
