import argparse
import bisect
import os
import subprocess
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path
import sys
import pandas as pd


def get_vcf_records_count(vcf_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if not (Path(f'{vcf_file}.tbi').exists() or Path(f'{vcf_file}.csi').exists()):
        os.system(f'tabix -f -p vcf {vcf_file}')
    process = subprocess.Popen(['bcftools', 'index', '-n', vcf_file],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        return 0
    return int(stdout)


def is_vcf_chrom_code_contains_chr(vcf_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if not (Path(f'{vcf_file}.tbi').exists() or Path(f'{vcf_file}.csi').exists()):
        os.system(f'tabix -f -p vcf {vcf_file}')
    process = subprocess.Popen(['tabix', '-l', vcf_file],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        return False
    return 'chr' in stdout.decode("UTF-8")


def find_closest_from_list(target_list, value):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    pos = bisect.bisect_left(target_list, value)
    if pos == 0:
        return target_list[0]
    if pos == len(target_list):
        return target_list[len(target_list) - 1]
    before = target_list[pos - 1]
    after = target_list[pos]
    if after - value < value - before:
        return after
    else:
        return before


def calc_for_chrom(full_gene_list_file, chrom, input_vcf, maf_thresh=0.1):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    start_time = datetime.now()
    chrom_code_in_vcf = f'chr{chrom}' if is_vcf_chrom_code_contains_chr(input_vcf) else chrom
    result_file = f'ld_r2_stats_{chrom_code_in_vcf}.tsv'
    print(f'Calculate r2 stats for chromosome {chrom} start at: {start_time}')
    gene_info_df = pd.read_table(full_gene_list_file,
                                 dtype={'start': 'Int64', 'end': 'Int64', 'position': 'Int64', 'chrom': 'string'})
    gene_info_df = gene_info_df[gene_info_df['chrom'] == str(chrom)]
    gene_info_df.drop(
        labels=gene_info_df[~gene_info_df['gene_type'].isin(['protein_coding', 'pseudogene', 'lincRNA'])].index,
        inplace=True)
    gene_info_df['chrom'] = gene_info_df['chrom'].astype(int)
    gene_info_df.sort_values(['chrom', 'start'], inplace=True)
    gene_info_df['center_snp'] = 'NA'
    gene_info_df['center_snp_pos'] = -1
    gene_info_df['r2_0_04_cnt'] = 0
    gene_info_df['r2_04_07_cnt'] = 0
    gene_info_df['r2_07_09_cnt'] = 0
    gene_info_df['r2_09_1_cnt'] = 0
    gene_len_thresh = 200000
    for idx, row in gene_info_df.iterrows():
        gene_len = row.loc['end'] - row.loc['start']
        distance_to_include = round((gene_len_thresh - gene_len) / 2) if gene_len < gene_len_thresh else 0
        start = max(row.loc['start'] - distance_to_include, 1)
        end = row.loc['end'] + distance_to_include
        gene_id = row.loc['gene_id']
        subset_vcf = f'subset_{chrom_code_in_vcf}.vcf.gz'
        with open(subset_vcf, mode='w') as f:
            proc1 = subprocess.Popen(['tabix', '-h', input_vcf, f'{chrom_code_in_vcf}:{start}-{end}'],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc2 = subprocess.Popen(['bcftools', 'view', '-m2', '-M2', '-v', 'snps', '-Ou'],
                                     stdin=proc1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc1.stdout.close()
            proc3 = subprocess.Popen(['bcftools', 'norm', '-d', 'any', '-Oz'],
                                     stdin=proc2.stdout, stdout=f, stderr=subprocess.PIPE)
            proc2.stdout.close()
            # bcftools output summary numbers to stderr,
            # so can not use stderr of proc3.communicate to check if subproces is success
            # https://github.com/samtools/bcftools/issues/850
            proc1.wait()
            proc2.wait()
            proc3.communicate()
        os.system(f'tabix -f -p vcf {subset_vcf}')
        if get_vcf_records_count(subset_vcf) == 0:
            print(f'warning: no variant in {gene_id}')
            continue
        # merge SNP rsid to target vcf
        print(f'merging SNP column to vcf data for gene {gene_id}')
        subset_snp_vcf = f'fill_snp_{chrom_code_in_vcf}.vcf'
        os.system(f'bcftools annotate -r {chrom_code_in_vcf} -I +"%CHROM:%POS\\_%REF\\_%FIRST_ALT" '
                  f'{subset_vcf} -Ov -o {subset_snp_vcf}')
        os.system(f'plink --silent --vcf {subset_vcf} --freq --out {chrom_code_in_vcf}_pre')
        maf_df = pd.read_table(f'{chrom_code_in_vcf}_pre.frq', sep=r'\s+')
        # include all variant that has maf>maf_thresh
        maf_df = maf_df[maf_df['MAF'] > maf_thresh]
        if maf_df.shape[0] < 2:
            print(f'warning: no variant with maf > {maf_thresh} for {gene_id}')
            continue
        vcf_df = pd.read_table(f'{subset_snp_vcf}', header=None, comment='#', usecols=[1, 2])
        vcf_df.columns = ['POS', 'SNP']
        maf_df = pd.merge(maf_df, vcf_df,
                          left_on='SNP',
                          right_on='SNP',
                          how='left')
        # 距离maf_df中心最近的且MAF>maf_thresh的点作为center_snp即gwas causal
        center_snp_pos = find_closest_from_list(maf_df['POS'],
                                                (maf_df.at[0, 'POS'] + maf_df.at[maf_df.shape[0] - 1, 'POS']) / 2)
        center_snp_pos_idx = vcf_df[vcf_df['POS'] == center_snp_pos].index
        center_snp = vcf_df.at[center_snp_pos_idx[0], 'SNP']
        os.system(f'plink --silent --vcf {subset_snp_vcf} '
                  f'--r2 inter-chr with-freqs --ld-window-r2 0 '
                  f'--ld-snp {center_snp} '
                  f'--out ld_{chrom_code_in_vcf}')
        ld_df = pd.read_table(f'ld_{chrom_code_in_vcf}.ld', sep=r'\s+', header=0)
        if ld_df.shape[0] == 0:
            print(f'warning: ld is empty in {gene_id}')
            continue
        # drop the row which represent r2 with self(r2=1)
        ld_df.drop(index=ld_df.loc[ld_df['SNP_B'] == ld_df['SNP_A']].index, inplace=True)
        ld_df.drop(index=ld_df.loc[ld_df['MAF_B'] <= maf_thresh].index, inplace=True)
        # keep variants within 50kb of gwas causal loci
        eqtl_gwas_causal_distance = 50000
        ld_df.drop(index=ld_df.loc[(ld_df['BP_A'] - ld_df['BP_B']).abs() > eqtl_gwas_causal_distance].index,
                   inplace=True)
        if ld_df.shape[0] == 0:
            print(f'warning: ld is empty after filtering in {gene_id}')
            continue
        gene_info_df.loc[idx, 'center_snp'] = center_snp
        gene_info_df.loc[idx, 'center_snp_pos'] = center_snp_pos
        gene_info_df.loc[idx, 'r2_0_04_cnt'] = ld_df['R2'][ld_df['R2'] <= 0.4].count()
        gene_info_df.loc[idx, 'r2_04_07_cnt'] = ld_df['R2'][(ld_df['R2'] > 0.4) & (ld_df['R2'] <= 0.7)].count()
        gene_info_df.loc[idx, 'r2_07_09_cnt'] = ld_df['R2'][(ld_df['R2'] > 0.7) & (ld_df['R2'] <= 0.9)].count()
        gene_info_df.loc[idx, 'r2_09_1_cnt'] = ld_df['R2'][ld_df['R2'] > 0.9].count()
    gene_info_df.to_csv(result_file, sep='\t', header=True, index=False)
    print(f'Calculate r2 stats for chromosome {chrom} completed, duration {datetime.now() - start_time}')


def run(src_vcf_dir, gene_list, start, end, maf_thresh, proc_count):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    result_file = 'ld_r2_stats.tsv'
    if proc_count > 22:
        proc_count = 22
    if proc_count < 1:
        proc_count = 1
    p = Pool(proc_count)
    for i in range(start, end + 1):
        p.apply_async(calc_for_chrom, args=(gene_list, i, os.path.join(src_vcf_dir, f'chr{i}.vcf.gz'), maf_thresh))
    p.close()
    p.join()
    result_list = []
    for f in os.listdir(os.getcwd()):
        if f.startswith('ld_r2_stats_'):
            result_list.append(pd.read_table(f, header=0, sep=r'\s+'))
    result_df = pd.concat(result_list)
    result_df.sort_values(['chrom', 'start'], inplace=True)
    result_df.to_csv(result_file, sep='\t', header=True, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--src_vcf_dir', dest='src_vcf_dir',
                        default='/Volumes/HD/biodata/colocalization-tools/raw/vcf/phased_hg38/EUR',
                        help='Directory of source vcf, vcf name must be chr{num}.vcf.gz')
    parser.add_argument('--gene_list', dest='gene_list',
                        default='/Volumes/HD/biodata/colocalization-tools/raw/genecode/gtf_all_gene.tsv',
                        help='Candidate gene list to calculate')
    parser.add_argument('--maf', dest='maf', default=0.1, type=float, help='MAF threshold to pick center SNP')
    parser.add_argument('--start_chr', dest='start_chr', default=1, type=int, choices=range(1, 23),
                        help='First chromosome to calculate')
    parser.add_argument('--end_chr', dest='end_chr', default=22, type=int, choices=range(1, 23),
                        help='Last chromosome to calculate, inclusive')
    parser.add_argument('--proc_cnt', dest='proc_cnt', default=4, type=int,
                        help='Parallel running process count')
    args = parser.parse_args()
    print(f'Accepted args:\n {args}')
    run(args.src_vcf_dir, args.gene_list, args.start_chr, args.end_chr, args.maf, args.proc_cnt)
