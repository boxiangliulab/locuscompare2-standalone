# Data preparation

You can put the prepared data anywhere and use soft links to the `data` directory

Here is an example:

If you have prepared all files in the `/locuscompare2file` directory of your computer, You can create a soft link by running the following command in the current folder:

```
ln -s /locuscompare2file/GWAS.tsv.gz
ln -s /locuscompare2file/eQTL.tsv.gz
ln -s /locuscompare2file/vcf_hg38/
ln -s /locuscompare2file/genecode.v26.basic.annotation.gtf.gz
ln -s /locuscompare2file/ld_block_loci_file
ln -s /locuscompare2file/TWAS_weights/
ln -s /locuscompare2file/predixcan_db/
```

