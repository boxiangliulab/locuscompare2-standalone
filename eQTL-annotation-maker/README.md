# eQTL annotation maker

Preparing probabilistic eQTL annotation; derive annotations based on your own eQTL data

## 1. Setup Environment

1) [dap](https://github.com/xqwen/dap)
2) [torus](https://github.com/xqwen/torus)
3) perl
4) [qtltools](https://github.com/qtltools/qtltools)
5) bgzip
6) gzip

## 2. Input file

+ [genotype file](http://vcftools.sourceforge.net/specs.html)

  Any genotype data needs to be specified using VCF or BCF file format(the third columns format should chr_22_17470194_G_C)

```
#CHROM	POS	        ID                      REF     ALT
22      17470194    chr_22_17470194_G_C         G       C
```

+ phenotype file

  Phenotype data are specified using an extended UCSC BED format. It is a standard BED file with some additional columns

```
#chr    start       end         pid             gid              strand  id2_0       id2_1       id2_2   
1       32936444    32964685    ENSG00000116514 ENSG00000116514     -    0.00669943  -0.200025   -0.491984
```

+ covariate file

  The COV file contains the covariate data in simple TXT format

```
SampleID    id2_0       id2_1
PC1         0.290556    0.240027
PC2         2.73356     0.605133
```

+ perl script files

  This will require the three perl file in the current working directory.
    1. [assemble.pl](https://github.com/yuer0608/dap/blob/master/gtex_v8_analysis/assemble.pl)
    2. [process.pl](https://github.com/yuer0608/dap/blob/master/gtex_v8_analysis/process.pl)
    3. [summarize_dap2enloc.pl](https://github.com/xqwen/fastenloc/blob/master/src/summarize_dap2enloc.pl)

## 3. Run

+ Run eQTL annotation maker in command line

```shell
python eqtl_annotation_maker --d path_dir --perl_dir perl_scripts_dir --g genotype_file --p phenotype_file --c covariate_file [--t tissue_name] 
```

## 4. Output

After script, a generate annotation vcf file in the working directory.