import os


def main(src_dir, dest_dir):
    template = '''input:
  eqtl:
    col_name_mapping:
      alt: alt
      beta: beta
      chrom: var_chrom
      gene_id: phe_id
      maf: maf
      position: var_from
      pvalue: nom_pval
      ref: ref
      se: se
      snp: var_id
    file: /Volumes/HD/biodata/sem/eqtl/sum/@.tsv.gz
    sample_size: 500
    tissue: test_tissue
    type: quant
  eqtl_finemapping_file: /Volumes/HD/biodata/sem/fastenloc/@.vcf.gz
  genecode: /Volumes/HD/biodata/colocalization-tools/raw/genecode/gencode.v26.basic.annotation.gtf.gz
  gwas:
    col_name_mapping:
      beta: beta
      chrom: var_chrom
      eaf: maf
      effect_allele: alt
      other_allele: ref
      position: var_from
      pvalue: nom_pval
      se: se
      snp: var_id
    file: /Volumes/HD/biodata/sem/gwas/sum/@.tsv.gz
    sample_size: 500
    trait: @
    type: quant
  ld_block_loci_file: /Volumes/HD/biodata/colocalization-tools/raw/loci/eur_hg38_ld_block.bed
  prediction_dir: /Volumes/HD/biodata/sem/pred_model
  twas_model_dir: /Volumes/HD/biodata/sem/twas_model
  vcf: /Volumes/HD/biodata/colocalization-tools/raw/vcf/phased_hg38
neighbour_snp_range: 25
p-value_threshold:
  eqtl: 0.01
  gwas: 1.0e-05
population: EUR
working_dir: /Volumes/HD/biodata/sem_out
'''
    for src_file in os.listdir(src_dir):
        if src_file.startswith('.') or os.path.isdir(src_file):
            continue
        if not src_file.endswith('.vcf.gz'):
            continue
        gene_id = src_file.strip('.vcf.gz')
        cfg = template.replace('@', gene_id)
        with open(os.path.join(dest_dir, f'{gene_id}.yml'), 'w') as out:
            out.write(cfg)


if __name__ == '__main__':
    main('/Volumes/HD/biodata/sem/gwas', '/Volumes/HD/gitrepo/configs/sem')
