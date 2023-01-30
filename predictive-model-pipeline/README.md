# TWAS Predictive Model Pipeline

Computing predictive models and putting them in a format that can be used to run [TWAS](http://gusevlab.org/projects/fusion/).
R script file FUSION.compute_weights.R comes from [TWAS github repo](https://github.com/gusevlab/fusion_twas), 

## 1. Setup Environment

1) Download [FUSION.compute_weights.R](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.compute_weights.R) and add its path to PATH env.
2) Required R packages: optparse, glmnet, methods and [plink2R](https://github.com/gabraham/plink2R).
3) Other softwares: [tabix](http://www.htslib.org/doc/tabix.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [plink](https://www.cog-genomics.org/plink/), [gcta_nr_robust](https://github.com/gusevlab/fusion_twas/blob/master/gcta_nr_robust) or original [gcta](https://yanglab.westlake.edu.cn/software/gcta).
4) If including bslmm/blup model, download and install [GEMMA](http://www.xzlab.org/software.html).
5) Python packages: pandas.

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
+ Note that since GEMMA requires results to go into an output subdirectory in working directory, the pipline will create a symbolic link to output subdirectory in the working directory no matter what the --output_dir param specifies.

```shell
python3 predictive-model-pipeline/pipeline.py --genotype genotpe_vcf_path --expression gene_exp_bed_path --tissue_name tissue1 --additional_compute_weight_params "--PATH_gcta path/to/gcta --PATH_plink path/to/plink" --output_dir output_directory
```

## 4. Output

Results pos file is placed in [output_di rectory], including pos file and *.RDat file
