#!/bin/sh

#variant filtering: MAF>1%; missingness <1%; HWE P-value > 1e-6
#https://github.com/gusevlab/fusion_twas/issues/28

for i in $(seq 1 22)
do
	plink --vcf chr${i}.vcf.gz --maf 0.01 --geno 0.01 --mind 0.01 --hwe 1e-06 --snps-only just-acgt --biallelic-only strict --make-bed --out chr${i}
	rm chr${i}.nosex
	rm chr${i}.log
done
