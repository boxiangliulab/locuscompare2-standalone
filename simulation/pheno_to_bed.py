import pandas as pd
import sys
from datetime import datetime
import random
import argparse
import os

def convert_gcta_pheno_to_bed(input_pheno_file, output_bed_file, chrom, start_pos=None, end_pos=None,
                              gene_name=None, strand=None):
    # Suppose that we simulate one gene on one chromosome per time
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    start_time = datetime.now()
    print(f'Converting phenotype file {input_pheno_file} to bed file {output_bed_file} start at: {start_time}')
    pheno_df = pd.read_table(input_pheno_file, sep=r'\s+', header=None)
    # Drop the first column, i.e. Family ID column
    pheno_df.drop(columns=0, inplace=True)
    # Reset column names after deleted FID column
    pheno_df.columns = range(len(pheno_df.columns))
    # Use IID column as index, IID will become header after transpose
    pheno_df.set_index(0, inplace=True)
    # Add prefix 6 columns
    result = pheno_df.T
    # the chr column should be exactly the same as CHROM field in genotype vcf
    result['chr'] = chrom
    result['start'] = 0 if start_pos is None else start_pos
    # (end - start) 必须 >= len, 注意bed的end是exclusive的
    # 所以需要大于等于len, 保证每条记录可以被分一个组, QTLTools --chunk的时候会用到
    result['end'] = pheno_df.shape[0] if end_pos is None else end_pos
    dummy_gene_id_suffix = random.randint(1, 99999)
    result['pid'] = f'gene{dummy_gene_id_suffix}' if gene_name is None else gene_name
    result['gid'] = f'gene{dummy_gene_id_suffix}' if gene_name is None else gene_name
    result['strand'] = random.choice(['+', '-']) if strand is None else strand
    prefix_col = ['chr', 'start', 'end', 'pid', 'gid', 'strand']
    result = result.reindex(columns=prefix_col + [col for col in result.columns if col not in prefix_col], copy=False)
    # QTLTools require a header row
    with open(output_bed_file, mode='w') as file_ptr:
        file_ptr.write(f'#')
    # Append content to output file using mode append, header will be prefixed with #
    result.to_csv(output_bed_file, sep='\t', header=True, index=False, mode='a')
    print(f'Process completed, duration {datetime.now() - start_time}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', dest='input_file',
                        required=True,
                        help='Input gcta simulation output .phen file')
    parser.add_argument('--output_file', dest='output_file',
                        required=True,
                        help='Output bed file path')
    parser.add_argument('--chromosome', dest='chromosome',
                        required=True,
                        help='Chromosome code, must be exactly the same as CHROM field in genotype vcf')
    parser.add_argument('--gene_name', dest='gene_name',
                        help='Gene name of the data to be converted')
    parser.add_argument('--strand', dest='strand', choices=['+', '-'],
                        help='Strand of the gene')
    parser.add_argument('--start_pos', dest='start_pos', type=int,
                        help='Start position of the gene')
    parser.add_argument('--end_pos', dest='end_pos', type=int,
                        help='End position of the gene')
    args = parser.parse_args()
    print(f'Accepted args:\n {args}')
    convert_gcta_pheno_to_bed(args.input_file, args.output_file, args.chromosome, args.gene_name, args.strand,
                              args.start_pos, args.end_pos)
