import argparse
import gzip
import logging
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

SNP_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


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


def parse_gtf(gtf_path, out_path):
    """
    Reference from parse_gtf.py
    The start column output is different from pyranges, start column will not minus 1 here
    """
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info('Parse gtf start')
    header_fields = ['chr', 'gene_id', 'gene_name', 'start', 'end', 'gene_type']
    gtf_compressed = gtf_path.endswith('.gz')
    with gzip.open(gtf_path, 'rb') if gtf_compressed else open(gtf_path) as gtf, open(out_path, 'w') as out:
        out_header = '\t'.join(header_fields) + '\n'
        out.write(out_header)
        for line in gtf:
            # Skip comments.
            if gtf_compressed:
                line = line.decode('UTF-8')
            if line[0] == '#':
                continue
            gene_fields = line.split('\t')
            # Exclude any rows not describing genes.
            if gene_fields[2] != 'gene':
                continue
            gene_attributes = gene_fields[-1].split('; ')
            attr_dict = dict(attribute.split() for attribute in gene_attributes if attribute)
            # Some gtf files may have chromosome number with 'chr' prefix.
            # We just need the number.
            chrom = gene_fields[0][3:] if "chr" in gene_fields[0] else gene_fields[0]
            # Extract rest of desired fields.
            start = gene_fields[3]
            end = gene_fields[4]
            gene_id = attr_dict['gene_id'].strip('"')
            gene_id = gene_id.split('.')[0]
            gene_name = attr_dict['gene_name'].strip('"')
            gene_type = attr_dict['gene_type'].strip('"')
            # Concatenate together and write out to file.
            out_line = '\t'.join([chrom, gene_id, gene_name, start, end, gene_type]) + '\n'
            out.write(out_line)
    logging.info('Parse gtf complete')


def create_snp_annot_from_vcf(input_vcf, output_dir, use_varid_in_covariances=True):
    """
    Reference from split_snp_annot_by_chr.py
    Output files are output_dir/snp_annot_chr{num}.tsv
    Note that chrom in varID must be prefixed with 'chr' since predixcan will add 'chr' prefix if varID ends with _b38,
    this will cause varID mismatch between model db and covariances, further causing predixcan reporting 0% snps used.
    see PredictionModel.Model.load_model
    """
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info('Create snp annotation from vcf start')
    if not os.path.exists(input_vcf) or os.path.getsize(input_vcf) <= 0:
        raise ValueError(f'{input_vcf} does not exist or is empty')
    header_fields = ['chromosome', 'pos', 'varID', 'ref_vcf', 'alt_vcf', 'rsid']
    vcf_df = pd.read_table(input_vcf, header=None, comment='#', usecols=[0, 1, 2, 3, 4],
                           dtype={0: 'category', 1: 'Int64', 2: 'string', 3: 'category', 4: 'category'})
    vcf_df.columns = ['chromosome', 'pos', 'rsid', 'ref_vcf', 'alt_vcf']
    # Drop non-autosome data
    vcf_df.drop(index=vcf_df[~vcf_df['chromosome'].isin([str(i) for i in range(1, 23)])].index,
                inplace=True)
    vcf_chr_notation = vcf_df['chromosome'].astype(str).str.contains('chr').any()
    vcf_df['varID'] = \
        vcf_df['chromosome'].astype(str) if vcf_chr_notation else 'chr' + vcf_df['chromosome'].astype(str) \
                                                                  + '_' + vcf_df['pos'].astype(str) \
                                                                  + '_' + vcf_df['ref_vcf'].astype(str) \
                                                                  + '_' + vcf_df['alt_vcf'].astype(str) \
                                                                  + '_b38'
    if use_varid_in_covariances:
        vcf_df['rsid'] = vcf_df['varID']
    vcf_df = vcf_df.reindex(columns=header_fields, copy=False)
    # Drop non-single letter polymorphisms
    indel_bool_series = (vcf_df['ref_vcf'].astype(str).str.len() != 1) | (vcf_df['alt_vcf'].astype(str).str.len() != 1)
    vcf_df.drop(labels=vcf_df[indel_bool_series].index, inplace=True)
    # Drop ambiguous strands data
    ambiguous_strand_series = (vcf_df['ref_vcf'] == vcf_df['alt_vcf'].map(SNP_COMPLEMENT))
    vcf_df.drop(labels=vcf_df[ambiguous_strand_series].index, inplace=True)
    # Drop missing rsid data
    vcf_df.drop(labels=vcf_df[vcf_df['rsid'] == '.'].index, inplace=True)
    if vcf_df.shape[0] == 0:
        raise ValueError(f'No eligible data in {input_vcf}, is ID column value missing(i.e. dot) in vcf?')
    grouped = vcf_df.groupby('chromosome')
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for name, group in grouped:
        chr_notation = f'{name}' if vcf_chr_notation else f'chr{name}'
        output_file = os.path.join(output_dir, f'snp_annot_{chr_notation}.tsv')
        group.to_csv(output_file, sep='\t', header=True, index=False, na_rep='NA')
    logging.info('Create snp annotation from vcf complete')


def create_genotype_from_vcf(input_vcf, output_dir):
    """
    Reference from split_genotype_by_chr.py
    Output files are output_dir/geno_chr{num}.tsv
    chrom in varID must be prefixed with 'chr'.
    """
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info('Create genotype files from vcf start')
    if not os.path.exists(input_vcf) or os.path.getsize(input_vcf) <= 0:
        raise ValueError(f'{input_vcf} does not exist or is empty')
    for line in gzip.open(input_vcf, 'rb'):
        line = line.decode('UTF-8').strip('\n')
        if line.startswith('#CHROM'):
            vcf_header = line.split('\t')
            break
    else:
        raise ValueError(f'{input_vcf} does not have a header row start with #CHROM')
    vcf_df = pd.read_table(input_vcf, header=None, comment='#')
    vcf_df.columns = vcf_header
    vcf_df.replace({'0/0': '0', '0/1': '1', '1/0': '1', '1/1': '2', './.': pd.NA}, inplace=True)
    vcf_chr_notation = vcf_df['#CHROM'].astype(str).str.contains('chr').any()
    # COUNTED column is major allele, .traw file
    vcf_df['varID'] = \
        vcf_df['#CHROM'].astype(str) if vcf_chr_notation else 'chr' + vcf_df['#CHROM'].astype(str) \
                                                              + '_' + vcf_df['POS'].astype(str) \
                                                              + '_' + vcf_df['REF'].astype(str) \
                                                              + '_' + vcf_df['ALT'].astype(str) \
                                                              + '_b38'
    indel_bool_series = (vcf_df['ALT'].str.len() != 1) | (vcf_df['REF'].str.len() != 1)
    vcf_df.drop(labels=vcf_df[indel_bool_series].index, inplace=True)
    ambiguous_strand_series = (vcf_df['ALT'] == vcf_df['REF'].map(SNP_COMPLEMENT))
    vcf_df.drop(labels=vcf_df[ambiguous_strand_series].index, inplace=True)
    vcf_df.drop_duplicates(subset='varID', inplace=True)
    if vcf_df.shape[0] == 0:
        raise ValueError(f'No eligible data in {input_vcf}')
    vcf_df.drop(columns=['POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
    grouped = vcf_df.groupby('#CHROM')
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for name, geno_group in grouped:
        geno_group.drop(columns='#CHROM', inplace=True)
        geno_group = geno_group.reindex(columns=['varID'] + [col for col in geno_group.columns if col not in ['varID']],
                                        copy=False)
        chr_notation = f'{name}' if vcf_chr_notation else f'chr{name}'
        output_file = os.path.join(output_dir, f'geno_{chr_notation}.tsv')
        geno_group.to_csv(output_file, sep='\t', header=True, index=False, na_rep='NA')
    logging.info('Create genotype files from vcf complete')


def preprocess_gene_exp(input_bed, output_for_peer, output_for_train,
                        gene_id_col_name=None, non_individual_col_names=None):
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    if non_individual_col_names is None:
        non_individual_col_names = ['#chr', 'start', 'end', 'pid', 'gid', 'strand']
    if gene_id_col_name is None:
        gene_id_col_name = 'pid'
    if gene_id_col_name in non_individual_col_names:
        non_individual_col_names.remove(gene_id_col_name)
    logging.info('Preprocess gene expression start')
    if not os.path.exists(input_bed) or os.path.getsize(input_bed) <= 0:
        raise ValueError(f'{input_bed} does not exist or is empty')
    bed_df = pd.read_table(input_bed, header=0)
    bed_df.drop(columns=non_individual_col_names, inplace=True)
    bed_df.set_index(gene_id_col_name, inplace=True)
    result = bed_df.T
    result.to_csv(output_for_peer, sep='\t', header=True, index=False)
    result.to_csv(output_for_train, sep='\t', header=True, index=True)
    logging.info('Preprocess gene expression complete')


def calc_peer_covariates(input_gene_exp, sample_size, individual_ids, output_dir, output_file_name):
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info('Calculate peer covariates start')
    if not os.path.exists(input_gene_exp) or os.path.getsize(input_gene_exp) <= 0:
        raise ValueError(f'{input_gene_exp} does not exist or is empty')
    if sample_size >= 350:
        n_factor = 60
    elif sample_size >= 250:
        n_factor = 45
    elif sample_size >= 150:
        n_factor = 30
    else:
        n_factor = 15
    os.system(f'peertool -f {input_gene_exp} -n {n_factor} --has_header -o {output_dir}')
    peer_factors = pd.read_table(os.path.join(output_dir, 'X.csv'), sep=',', header=None)
    peer_factors.columns = individual_ids
    peer_factors.to_csv(os.path.join(output_dir, output_file_name), sep='\t', header=True, index=True)
    # residuals are skipped
    logging.info('Calculate peer covariates complete')


def train(expression_file, gene_annot_file, covariates_file, geno_snp_input_dir,
          summary_dir, weights_dir, covariances_dir, sim_data=False):
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info('Training model start')
    process_cnt = os.cpu_count()
    if process_cnt > 22:
        process_cnt = 22
    if process_cnt < 1:
        process_cnt = 1
    p = Pool(process_cnt)
    for chrom in range(1, 23):
        p.apply_async(train_for_chrom, args=(chrom, expression_file, gene_annot_file, covariates_file,
                                             geno_snp_input_dir, summary_dir, weights_dir, covariances_dir, sim_data))
    p.close()
    p.join()
    cov_df_list = []
    result_cov_file = 'covariances.txt.gz'
    for cov_file in os.listdir(covariances_dir):
        if cov_file.endswith('_covariances.txt'):
            cov_df_list.append(pd.read_table(os.path.join(covariances_dir, cov_file), sep=' ', header=0))
    cov_df = pd.concat(cov_df_list)
    cov_df.to_csv(os.path.join(covariances_dir, result_cov_file), sep=' ', header=True, index=False)
    logging.info(f'Training model complete')


def train_for_chrom(chrom, expression_file, gene_annot_file, covariates_file, geno_snp_input_dir,
                    summary_dir, weights_dir, covariances_dir, sim_data=False):
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info(F'Training model for chromosome {chrom}')
    rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'train_elnet.R')
    snp_annot = os.path.join(geno_snp_input_dir, f'snp_annot_chr{chrom}.tsv')
    genotype = os.path.join(geno_snp_input_dir, f'geno_chr{chrom}.tsv')
    if not os.path.exists(snp_annot) or os.path.getsize(snp_annot) <= 0:
        logging.info(f'snp_annot for chrom {chrom} does exist, skipping')
        return
    if not os.path.exists(genotype) or os.path.getsize(genotype) <= 0:
        logging.info(f'genotype for chrom {chrom} does exist, skipping')
        return
    os.system(f'Rscript --no-save --no-restore {rscript_path} {snp_annot} {gene_annot_file}  {genotype} '
              f'{expression_file} {covariates_file} {chrom} {summary_dir} {weights_dir} {covariances_dir} {sim_data}')


def generate_db(summary_dir, weights_dir, output_db):
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info('Generate db start')
    rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'gen_db.R')
    os.system(f'Rscript --no-save --no-restore {rscript_path} {summary_dir} {weights_dir} {output_db}')
    logging.info('Generate db complete')


def run(input_geno_vcf=None, use_varid_in_covariances=True, input_exp_bed=None, exp_gene_id_col_name=None,
        exp_non_individual_col_names=None, input_gene_code=None, sim_data=False, output_dir=None):
    __init_logger()
    print(os.path.basename(__file__))
    print(sys._getframe().f_code.co_name)
    logging.info('Training pipline start')
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    genecode = os.path.join(output_dir, 'genecode_parsed.tsv')
    parse_gtf(input_gene_code, genecode)
    create_snp_annot_from_vcf(input_geno_vcf, output_dir, use_varid_in_covariances)
    create_genotype_from_vcf(input_geno_vcf, output_dir)
    # NOTE: peer input file extension: .csv(comma separated) or .tab(tab separated),
    # file name with other extension will be assumed separated by single space
    input_for_peer = os.path.join(output_dir, 'peer_input.tab')
    expression_file = os.path.join(output_dir, 'expression_input.tsv')
    preprocess_gene_exp(input_exp_bed, input_for_peer, expression_file, exp_gene_id_col_name,
                        exp_non_individual_col_names)
    # NOTE that input_for_train has row indexes
    peer_cov_file_name = 'covariates.tsv'
    individual_ids = pd.read_table(expression_file, header=None, skiprows=1, usecols=[0])
    calc_peer_covariates(input_for_peer, individual_ids.shape[0], individual_ids[0], output_dir, peer_cov_file_name)
    peer_covariates = os.path.join(output_dir, peer_cov_file_name)
    summary_dir = os.path.join(output_dir, 'summary')
    weights_dir = os.path.join(output_dir, 'weights')
    covariances_dir = os.path.join(output_dir, 'covariances')
    train(expression_file, genecode, peer_covariates, output_dir, summary_dir, weights_dir, covariances_dir, sim_data)
    output_db = os.path.join(output_dir, 'db', 'result.db')
    generate_db(summary_dir, weights_dir, output_db)
    logging.info('Training pipline complete')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genotype', dest='genotype_vcf',
                        help='eQTL genotype vcf format')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--use_varid_in_covariances', dest='use_varid_in_covariances',
                       action="store_true",
                       help='Use varID in result covariances file instead of rsid')
    group.add_argument('--use_rsid_in_covariances', dest='use_rsid_in_covariances',
                       action="store_true",
                       help='Use rsid in result covariances file instead of varID')
    parser.add_argument('--expression', dest='gene_expression',
                        help='Gene expression in bed format')
    parser.add_argument('--exp_gene_id_col_name', dest='exp_gene_id_col_name',
                        help='Gene id column name in Gene expression bed format file')
    parser.add_argument('--exp_non_individual_col_names', dest='exp_non_individual_col_names', type=str,
                        nargs='*',
                        help='Column names in Gene expression bed format file, '
                             'except individual column and gene id column')
    parser.add_argument('--genecode', dest='genecode',
                        help='Genecode gtf')
    parser.add_argument('--output_dir', dest='output_dir', help='Output directory')
    parser.add_argument('--sim_data', dest='sim_data',
                        action="store_true", help='data is simulation data')
    args = parser.parse_args()
    print(f'Accepted args:\n {args}')
    use_varid = args.use_varid_in_covariances if args.use_varid_in_covariances else not args.use_rsid_in_covariances
    run(input_geno_vcf=args.genotype_vcf, use_varid_in_covariances=use_varid,
        input_exp_bed=args.gene_expression, exp_gene_id_col_name=args.exp_gene_id_col_name,
        exp_non_individual_col_names=args.exp_non_individual_col_names, input_gene_code=args.genecode,
        sim_data=args.sim_data, output_dir=args.output_dir)
