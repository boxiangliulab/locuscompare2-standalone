import argparse
import logging
import os
import shutil
import sys
import uuid
from pathlib import Path

import pandas as pd


def __init_logger():
    log_format = '%(asctime)s::%(levelname)s::%(name)s::%(filename)s::%(lineno)d::%(message)s'
    stdout_hd = logging.StreamHandler(sys.stdout)
    stdout_hd.setLevel(logging.INFO)
    logging.root.handlers = []
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            stdout_hd,
        ]
    )


def run(rscript_path, input_geno_vcf=None, input_exp_bed=None, exp_gene_id_col_name=None, chr_col_name=None,
        start_col_name=None, end_col_name=None, strand_col_name=None, exp_non_individual_col_names=None,
        tissue_name=None, additional_compute_weight_params='', cis_len=1000, output_dir=None):
    __init_logger()
    if exp_non_individual_col_names is None:
        exp_non_individual_col_names = ['#chr', 'start', 'end', 'pid', 'gid', 'strand']
    if exp_gene_id_col_name is None:
        exp_gene_id_col_name = 'pid'
    if chr_col_name is None:
        chr_col_name = '#chr'
    if start_col_name is None:
        start_col_name = 'start'
    if end_col_name is None:
        end_col_name = 'end'
    if strand_col_name is None:
        strand_col_name = 'strand'
    if cis_len <= 0:
        cis_len = 1000
    if rscript_path is None:
        rscript_path = shutil.which('FUSION.compute_weights.R')
    logging.info('Computing pipline start')
    tissue_rdat_dir = os.path.join(output_dir, tissue_name)
    if rscript_path is None or not Path(rscript_path).exists():
        raise ValueError(f'Full path to FUSION.compute_weights.R is required or add it to PATH')
    Path(tissue_rdat_dir).mkdir(parents=True, exist_ok=True)
    pos_file = os.path.join(output_dir, f'{tissue_name}.pos')
    if Path(pos_file).exists():
        os.remove(pos_file)
    if not os.path.exists(input_exp_bed) or os.path.getsize(input_exp_bed) <= 0:
        raise ValueError(f'{input_exp_bed} does not exist or is empty')
    if not os.path.exists(input_geno_vcf) or os.path.getsize(input_geno_vcf) <= 0:
        raise ValueError(f'{input_geno_vcf} does not exist or is empty')
    bed_df = pd.read_table(input_exp_bed, header=0)
    if exp_gene_id_col_name not in exp_non_individual_col_names:
        exp_non_individual_col_names.append(exp_gene_id_col_name)
    if chr_col_name not in exp_non_individual_col_names:
        exp_non_individual_col_names.append(chr_col_name)
    if start_col_name not in exp_non_individual_col_names:
        exp_non_individual_col_names.append(start_col_name)
    if end_col_name not in exp_non_individual_col_names:
        exp_non_individual_col_names.append(end_col_name)
    if strand_col_name not in exp_non_individual_col_names:
        exp_non_individual_col_names.append(strand_col_name)
    individual_col_names = [col for col in bed_df.columns if col not in exp_non_individual_col_names]
    #  Generate a symbolic link to the gemma output dir,
    #  this is a workaround because GEMMA requires results to go into an output subdirectory in working directory.
    if not Path('output').exists():
        os.system(f'ln -s ./ output')
    with open(pos_file, 'w') as pos_f:
        pos_f.write(f'PANEL\tWGT\tID\tCHR\tP0\tP1\tN\n')
        for idx, row in bed_df.iterrows():
            gene_id = row.loc[exp_gene_id_col_name]
            chrom = row.loc['#chr']
            start = row.loc['start']
            end = row.loc['end']
            tss = start if row.loc['strand'] == '+' else end
            geno_start = max(0, tss - cis_len * 1000)
            geno_end = tss + cis_len * 1000
            if not os.path.exists(f'{input_geno_vcf}.tbi'):
                os.system(f'tabix -f -p vcf {input_geno_vcf}')
            subset_vcf = os.path.join(output_dir, f'{gene_id}_geno.vcf')
            os.system(f'tabix -h {input_geno_vcf} {chrom}:{geno_start}-{geno_end} > {subset_vcf}'
                      f'&& bgzip -f {subset_vcf}')
            subset_genotype_file = f'{subset_vcf}.gz'
            if not os.path.exists(subset_genotype_file) or os.path.getsize(subset_genotype_file) <= 0:
                logging.warning(f'Extracting genotypes for the cis locus of gene {gene_id} failed, skipping')
                continue
            phenotype_file = os.path.join(output_dir, f'{gene_id}_phen.tsv')
            phenotype_df = bed_df.loc[[idx], individual_col_names].T
            phenotype_df['fid'] = phenotype_df.index
            phenotype_df['iid'] = phenotype_df.index
            phenotype_df = phenotype_df.reindex(
                columns=['fid', 'iid'] + [col for col in phenotype_df.columns if col not in ['fid', 'iid']], copy=False)
            phenotype_df.to_csv(phenotype_file, sep='\t', index=False, header=False)
            phen_bin_prefix = os.path.join(output_dir, f'{gene_id}_phen_bin')
            os.system(f'plink --vcf {subset_genotype_file} --double-id '
                      f'--pheno {phenotype_file} --make-bed --out {phen_bin_prefix}')
            phen_bed = f'{phen_bin_prefix}.bed'
            if not os.path.exists(phen_bed) or os.path.getsize(phen_bed) <= 0:
                logging.warning(f'Plink binary phenotype file of gene {gene_id} generation failed, skipping')
                os.system(f'rm {subset_genotype_file}')
                os.system(f'rm {phenotype_file}')
                continue
            # tmp file must be in working directory and can only be relative path
            # since GEMMA will write output to ./output/{tmp_file}.X.Y, else GEMMA will fail to write output
            tmp_file = str(uuid.uuid4())
            weight_out_prefix = os.path.join(tissue_rdat_dir, f'{gene_id}')
            if additional_compute_weight_params is None:
                additional_compute_weight_params = ''
            os.system(f'Rscript --no-save --no-restore {rscript_path} '
                      f'--save_hsq --PATH_gcta gcta --bfile {phen_bin_prefix} '
                      f'--tmp {tmp_file} '
                      f'{additional_compute_weight_params} '
                      f'--out {weight_out_prefix}')
            weight_abs_path = f'{weight_out_prefix}.wgt.RDat'
            if not os.path.exists(weight_abs_path) or os.path.getsize(weight_abs_path) <= 0:
                logging.warning(f'No weight data generated for gene {gene_id}, skipping')
                os.system(f'rm {subset_genotype_file}')
                os.system(f'rm {phen_bin_prefix}.*')
                os.system(f'rm {phenotype_file}')
                continue
            weight_relative_path = os.path.join(tissue_name, f'{gene_id}.wgt.RDat')
            pos_f.write(f'{tissue_name}\t{weight_relative_path}\t'
                        f'{gene_id}\t{chrom}\t{start}\t{end}\t{len(individual_col_names)}\n')
            os.system(f'rm {subset_genotype_file}')
            os.system(f'rm {phen_bin_prefix}.*')
            os.system(f'rm {phenotype_file}')
    logging.info('Computing pipline complete')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rscript_path', dest='rscript_path',
                        help='Full path to FUSION.compute_weights.R if it is not added to PATH')
    parser.add_argument('--genotype', dest='genotype_vcf',
                        help='eQTL genotype vcf format, required')
    parser.add_argument('--expression', dest='gene_expression',
                        help='Gene expression in bed format, required')
    parser.add_argument('--exp_gene_id_col_name', dest='exp_gene_id_col_name',
                        help='Gene id column name in Gene expression bed format file')
    parser.add_argument('--exp_chr_col_name', dest='exp_chr_col_name',
                        help='Chromosome column name in Gene expression bed format file')
    parser.add_argument('--exp_start_col_name', dest='exp_start_col_name',
                        help='Chromosome start position column name in Gene expression bed format file')
    parser.add_argument('--exp_end_col_name', dest='exp_end_col_name',
                        help='Chromosome end position column name in Gene expression bed format file')
    parser.add_argument('--exp_strand_col_name', dest='exp_strand_col_name',
                        help='Gene strand column name in Gene expression bed format file')
    parser.add_argument('--exp_non_individual_col_names', dest='exp_non_individual_col_names', type=str,
                        nargs='*',
                        help='Column names in Gene expression bed format file, '
                             'except individual column and gene id column')
    parser.add_argument('--tissue_name', dest='tissue_name',
                        help='Tissue name, required')
    parser.add_argument('--additional_compute_weight_params', dest='additional_compute_weight_params',
                        help='Additional params for FUSION.compute_weights.R quoted by double quotes '
                             'e.g.: "--PATH_plink path/to/plink --PATH_gemma path/to/gemma"')
    parser.add_argument('--cis_len', dest='cis_len', default=1000, type=int,
                        help='cis locus length of a gene, unit kb, default 1000, '
                             'i.e. how long to extend on both sides of a gene to extract genotype data')
    parser.add_argument('--output_dir', dest='output_dir', help='Output directory, required')
    args = parser.parse_args()
    print(f'Accepted args:\n {args}')
    run(rscript_path=args.rscript_path, input_geno_vcf=args.genotype_vcf, input_exp_bed=args.gene_expression,
        exp_gene_id_col_name=args.exp_gene_id_col_name, chr_col_name=args.exp_chr_col_name,
        start_col_name=args.exp_start_col_name, end_col_name=args.exp_end_col_name,
        strand_col_name=args.exp_strand_col_name, exp_non_individual_col_names=args.exp_non_individual_col_names,
        tissue_name=args.tissue_name, additional_compute_weight_params=args.additional_compute_weight_params,
        cis_len=args.cis_len, output_dir=args.output_dir)
