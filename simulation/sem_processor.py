import argparse
import os
from datetime import datetime

import pandas as pd


def convert_sbams_to_merged(sbams_dir, dest_dir, prefix='out',
                            genecode='/ebs2/linyi/00_data/locuscompare2/raw/genecode/v26/gtf_all_gene.tsv',
                            individual_num=500,
                            dbsnp_ref_path='/ebs2/linyi/00_data/locuscompare2/raw/genotype/phased_hg38/EUR/chr5.vcf.gz'
                            ):
    pre_vcf_out = os.path.join(dest_dir, f'pre_{prefix}.vcf')
    bed_out = os.path.join(dest_dir, f'{prefix}.bed')
    ind_ids = [f'id{k}' for k in range(1, individual_num + 1)]
    genecode_df = pd.read_table(genecode, usecols=['start', 'end', 'strand', 'gene_id'])
    genecode_df.rename({'gene_id': 'gid'}, axis='columns', inplace=True)
    iid_mapper = {k: f'id{k - 2}' for k in range(3, individual_num + 3)}
    phen_df_list = []
    dosage_df_list = []
    for sbam in os.listdir(sbams_dir):
        if not sbam.endswith('.sbams.dat'):
            continue
        sbam_path = os.path.join(sbams_dir, sbam)

        dosage_df = pd.read_table(sbam_path, sep=r'\s+', header=None, skiprows=1,
                                  dtype={k: 'str' for k in range(3, individual_num + 3)})

        split_df = dosage_df[1].str.split("_", n=4, expand=True)
        chrom = int(split_df[0].loc[0].strip('chr'))
        dosage_df.drop(columns=[0, 2], inplace=True)
        dosage_df['FORMAT'] = 'GT'
        dosage_df['INFO'] = '.'
        dosage_df['FILTER'] = '.'
        dosage_df['QUAL'] = '.'
        dosage_df['ALT'] = split_df[3]
        dosage_df['REF'] = split_df[2]
        dosage_df['POS'] = split_df[1].astype(int)
        del split_df
        dosage_df['#CHROM'] = chrom
        dosage_mapper = iid_mapper.copy()
        dosage_mapper[1] = 'ID'
        dosage_df.rename(dosage_mapper, axis='columns', inplace=True)
        dosage_df = dosage_df.reindex(
            columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + ind_ids, copy=False)
        dosage_df.replace({'0': '0/0', '1': '0/1', '2': '1/1'}, inplace=True)
        dosage_df.fillna('./.', inplace=True)
        dosage_df_list.append(dosage_df)

        phen_df = pd.read_table(sbam_path, sep=r'\s+', header=None, nrows=1)
        # gene_id = phen_df.iat[0, 1]
        phen_df.drop(columns=[0, 2], inplace=True)
        phen_mapper = iid_mapper.copy()
        phen_mapper[1] = 'gid'
        phen_df.rename(phen_mapper, axis='columns', inplace=True)
        phen_df = pd.merge(left=phen_df, right=genecode_df, left_on='gid', right_on='gid', how='left')
        phen_df['#chr'] = chrom
        phen_df['pid'] = phen_df['gid']
        phen_df = phen_df.reindex(columns=['#chr', 'start', 'end', 'pid', 'gid', 'strand'] + ind_ids, copy=False)
        phen_df_list.append(phen_df)
    with open(pre_vcf_out, mode='w') as vcf_file:
        vcf_file.write(f'##fileformat=VCFv4.2\n')
        vcf_file.write(f'##FILTER=<ID=PASS,Description="All filters passed">\n')
        vcf_file.write(f'##fileDate={datetime.now().strftime("%Y%m%d")}\n')
        vcf_file.write(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf_file.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ind_ids) + '\n')
    with open(bed_out, mode='w') as bed_file:
        bed_file.write(f'#chr\tstart\tend\tpid\tgid\tstrand\t' + '\t'.join(ind_ids) + '\n')
    phen = pd.concat(phen_df_list)
    phen.drop_duplicates(subset='gid', inplace=True)
    phen.sort_values('start', inplace=True)
    phen.to_csv(bed_out, sep='\t', header=False, index=False, mode='a')
    del phen
    dosage = pd.concat(dosage_df_list)
    dosage.drop_duplicates(subset=['#CHROM', 'POS'], inplace=True)
    dosage.sort_values('POS', inplace=True)
    dosage.to_csv(pre_vcf_out, sep='\t', header=False, index=False, mode='a')
    del dosage
    #
    vcf_out = os.path.join(dest_dir, f'{prefix}.vcf.gz')
    os.system(f'bgzip -f {pre_vcf_out} && tabix -f -p vcf {pre_vcf_out}.gz && '
              f'bcftools annotate -a {dbsnp_ref_path} -c ID {pre_vcf_out}.gz -Oz -o {vcf_out} && tabix -f -p vcf {vcf_out}')
    os.remove(f'{pre_vcf_out}.gz')
    os.remove(f'{pre_vcf_out}.gz.tbi')

    os.system(f'bgzip -f {bed_out} && tabix -f -p bed {bed_out}.gz')
    pca_out = f'{prefix}.pca'
    os.system(
        f'QTLtools pca --bed {bed_out}.gz --out {prefix}_raw_pca1 --center --scale && head -4 {prefix}_raw_pca1 > {pca_out} && gzip -k {pca_out}')
    os.remove(f'{prefix}_raw_pca1')
    if not os.path.exists(pca_out):
        raise ValueError('pca not generated')
    pre_sum_out = f'{prefix}_sum'
    os.system(
        f'QTLtools cis --vcf {vcf_out} --bed {bed_out}.gz --cov {pca_out} --nominal 1 --normal --std-err --out {pre_sum_out}')
    if not os.path.exists(pre_sum_out):
        raise ValueError('sumstats not generated')
    sumstats_df = pd.read_table(pre_sum_out, header=None, sep=r'\s+')
    if sumstats_df.shape[0] == 0:
        raise ValueError('sumstats is empty')
    os.remove(pre_sum_out)
    sumstats_df.columns = ['phe_id', 'phe_chrom', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var',
                           'var_id', 'var_chrom', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'beta', 'se',
                           'best_hit']
    os.system(f'plink --vcf {vcf_out} --freq --double-id --out  {prefix}_maf')
    maf_df = pd.read_table(f'{prefix}_maf.frq', sep=r'\s+', usecols=['SNP', 'MAF'])
    os.remove(f'{prefix}_maf.frq')
    os.remove(f'{prefix}_maf.nosex')
    os.remove(f'{prefix}_maf.log')
    maf_df.columns = ['var_id', 'maf']
    sumstats_df = pd.merge(sumstats_df, maf_df,
                           left_on='var_id',
                           right_on='var_id',
                           how='left')
    # merge REF, ALT from genotype vcf
    genotype_df = pd.read_table(vcf_out, header=None, comment='#', usecols=[2, 3, 4])
    genotype_df.columns = ['var_id', 'ref', 'alt']
    sumstats_df = pd.merge(sumstats_df, genotype_df,
                           left_on='var_id',
                           right_on='var_id',
                           how='left')
    sum_out = os.path.join(dest_dir, f'{prefix}.tsv.gz')
    sumstats_df.to_csv(sum_out, sep='\t', header=True, index=False)


def convert_sbams_to_split(sbams_dir, dest_dir,
                           genecode='/ebs2/linyi/00_data/locuscompare2/raw/genecode/v26/gtf_all_gene.tsv',
                           individual_num=500,
                           dbsnp_ref_path='/ebs2/linyi/00_data/locuscompare2/raw/genotype/phased_hg38/EUR/chr5.vcf.gz'):
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    sum_dir = os.path.join(dest_dir, 'sum')
    if not os.path.exists(sum_dir):
        os.makedirs(sum_dir)
    ind_ids = [f'id{k}' for k in range(1, individual_num + 1)]
    genecode_df = pd.read_table(genecode, usecols=['start', 'end', 'strand', 'gene_id'])
    genecode_df.rename({'gene_id': 'gid'}, axis='columns', inplace=True)
    iid_mapper = {k: f'id{k - 2}' for k in range(3, individual_num + 3)}
    for sbam in os.listdir(sbams_dir):
        if not sbam.endswith('.sbams.dat'):
            continue
        sbam_path = os.path.join(sbams_dir, sbam)

        dosage_df = pd.read_table(sbam_path, sep=r'\s+', header=None, skiprows=1,
                                  dtype={k: 'str' for k in range(3, individual_num + 3)})

        split_df = dosage_df[1].str.split("_", n=4, expand=True)
        chrom = int(split_df[0].loc[0].strip('chr'))
        dosage_df.drop(columns=[0, 2], inplace=True)
        dosage_df['FORMAT'] = 'GT'
        dosage_df['INFO'] = '.'
        dosage_df['FILTER'] = '.'
        dosage_df['QUAL'] = '.'
        dosage_df['ALT'] = split_df[3]
        dosage_df['REF'] = split_df[2]
        dosage_df['POS'] = split_df[1].astype(int)
        del split_df
        dosage_df['#CHROM'] = chrom
        dosage_mapper = iid_mapper.copy()
        dosage_mapper[1] = 'ID'
        dosage_df.rename(dosage_mapper, axis='columns', inplace=True)
        dosage_df = dosage_df.reindex(
            columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + ind_ids, copy=False)
        dosage_df.replace({'0': '0/0', '1': '0/1', '2': '1/1'}, inplace=True)
        dosage_df.fillna('./.', inplace=True)

        phen_df = pd.read_table(sbam_path, sep=r'\s+', header=None, nrows=1)
        gene_id = phen_df.iat[0, 1]
        pre_vcf_out = os.path.join(dest_dir, f'pre_{gene_id}.vcf')
        vcf_out = os.path.join(dest_dir, f'{gene_id}.vcf.gz')
        with open(pre_vcf_out, mode='w') as vcf_file:
            vcf_file.write(f'##fileformat=VCFv4.2\n')
            vcf_file.write(f'##FILTER=<ID=PASS,Description="All filters passed">\n')
            vcf_file.write(f'##fileDate={datetime.now().strftime("%Y%m%d")}\n')
            vcf_file.write(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            vcf_file.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ind_ids) + '\n')
        dosage_df.drop_duplicates(subset=['#CHROM', 'POS'], inplace=True)
        dosage_df.sort_values('POS', inplace=True)
        min_pos = dosage_df['POS'].min()
        max_pos = dosage_df['POS'].max()
        dosage_df.to_csv(pre_vcf_out, sep='\t', header=False, index=False, mode='a')
        del dosage_df
        os.system(f'bgzip -f {pre_vcf_out} && tabix -f -p vcf {pre_vcf_out}.gz && '
                  f'bcftools annotate -r {chrom}:{min_pos}-{max_pos} -a {dbsnp_ref_path} -c ID {pre_vcf_out}.gz -Oz -o {vcf_out} && tabix -f -p vcf {vcf_out}')
        os.remove(f'{pre_vcf_out}.gz')
        os.remove(f'{pre_vcf_out}.gz.tbi')
        phen_df.drop(columns=[0, 2], inplace=True)
        phen_mapper = iid_mapper.copy()
        phen_mapper[1] = 'gid'
        phen_df.rename(phen_mapper, axis='columns', inplace=True)
        phen_df = pd.merge(left=phen_df, right=genecode_df, left_on='gid', right_on='gid', how='left')
        phen_df['#chr'] = chrom
        phen_df['pid'] = phen_df['gid']
        phen_df = phen_df.reindex(columns=['#chr', 'start', 'end', 'pid', 'gid', 'strand'] + ind_ids, copy=False)
        bed_out = os.path.join(dest_dir, f'{gene_id}.bed')
        with open(bed_out, mode='w') as bed_file:
            bed_file.write(f'#chr\tstart\tend\tpid\tgid\tstrand\t' + '\t'.join(ind_ids) + '\n')
        phen_df.drop_duplicates(subset='gid', inplace=True)
        phen_df.sort_values('start', inplace=True)
        phen_df.to_csv(bed_out, sep='\t', header=False, index=False, mode='a')
        del phen_df
        os.system(f'bgzip -f {bed_out} && tabix -f -p bed {bed_out}.gz')
        pre_sum_out = f'{gene_id}_sum'
        os.system(
            f'QTLtools cis --vcf {vcf_out} --bed {bed_out}.gz --nominal 1 --normal --std-err --out {pre_sum_out}')
        if not os.path.exists(pre_sum_out):
            print(f'Warning: QTL analysis result of {gene_id} does not exist')
            continue
        ana_result_df = pd.read_table(pre_sum_out, header=None, sep=r'\s+')
        os.remove(pre_sum_out)
        if ana_result_df.shape[0] == 0:
            print(f'Warning: QTL analysis result of {gene_id} is empty')
            continue
        ana_result_df.columns = ['phe_id', 'phe_chrom', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis',
                                 'dist_phe_var',
                                 'var_id', 'var_chrom', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'beta', 'se',
                                 'best_hit']
        os.system(f'plink --vcf {vcf_out} --freq --double-id --out  {gene_id}_maf')
        eqtl_maf_df = pd.read_table(f'{gene_id}_maf.frq', sep=r'\s+', usecols=['SNP', 'MAF'])
        os.remove(f'{gene_id}_maf.frq')
        os.remove(f'{gene_id}_maf.nosex')
        os.remove(f'{gene_id}_maf.log')
        eqtl_maf_df.columns = ['var_id', 'maf']
        ana_result_df = pd.merge(ana_result_df, eqtl_maf_df,
                                 left_on='var_id',
                                 right_on='var_id',
                                 how='left')
        # merge REF, ALT from genotype vcf
        genotype_df = pd.read_table(vcf_out, header=None, comment='#', usecols=[2, 3, 4])
        genotype_df.columns = ['var_id', 'ref', 'alt']
        ana_result_df = pd.merge(ana_result_df, genotype_df,
                                 left_on='var_id',
                                 right_on='var_id',
                                 how='left')
        sum_out = os.path.join(sum_dir, f'{gene_id}.tsv.gz')
        ana_result_df.to_csv(sum_out, sep='\t', header=True, index=False)


def split_eqtl_annotation(eqtl_annot='/Volumes/HD/biodata/sem/sim_scheme/sim_rst/coloc/eqtl.vcf.gz',
                          dest_dir='/Volumes/HD/biodata/sem/fastenloc/',
                          dbsnp_ref_path='/Volumes/HD/biodata/colocalization-tools/raw/vcf/phased_hg38/EUR/chr5.vcf.gz'):
    df = pd.read_table(eqtl_annot, header=None)
    split_df = df[2].str.split('_', expand=True, n=1)
    df['gene_id'] = split_df[0]
    df[2] = df[2].str.slice(16, -4)
    df[0] = df[0].str.strip('chr')
    df.columns = ['chrom', 'pos', 'var_id', 'ref', 'alt', 'weight', 'gene_id']
    cols = df.columns.to_list()
    cols.remove('gene_id')
    vcf_df = pd.read_table(dbsnp_ref_path, header=None, comment='#', usecols=[0, 1, 2, 3, 4],
                           dtype={0: 'category', 1: 'Int64'})
    vcf_df.columns = ['chrom', 'pos', 'rsid', 'ref', 'alt']
    vcf_df['var_id'] = 'chr' + vcf_df['chrom'].astype(str) + '_' + \
                       vcf_df['pos'].astype(str) + '_' + \
                       vcf_df['ref'].astype(str) + '_' + vcf_df['alt'].astype(str)
    vcf_df.drop(columns=['chrom', 'pos', 'ref', 'alt'], inplace=True)
    merged_df = pd.merge(left=df, right=vcf_df, left_on='var_id', right_on='var_id', how='left')
    merged_df['var_id'] = merged_df['rsid']
    for name, group in merged_df.groupby('gene_id'):
        group.replace({'@=': '@test_tissue='}, inplace=True, regex=True)
        group[cols].to_csv(os.path.join(dest_dir, f'{name}.vcf.gz'), sep='\t', header=False,
                           index=False)


def update_eqtl_annotation(eqtl_annot='/Volumes/HD/biodata/sem/sim_scheme/sim_rst/coloc/eqtl.vcf.gz',
                           dest_file='/Volumes/HD/biodata/sem/fastenloc/sem_eqtl_annot.vcf.gz'):
    df = pd.read_table(eqtl_annot, header=None)
    df[0] = df[0].str.strip('chr')
    df.columns = ['chrom', 'pos', 'var_id', 'ref', 'alt', 'weight']
    df.replace({'@=': '@test_tissue='}, inplace=True, regex=True)
    df.to_csv(dest_file, sep='\t', header=False, index=False)


def merge_sum(srd_dir='/Volumes/HD/biodata/sem/gwas/sum',
              dest_file='/Volumes/HD/biodata/sem/merged_gwas_sum_stats.tsv.gz'):
    df_list = []
    for gene_sum in os.listdir(srd_dir):
        if not gene_sum.startswith('ENSG'):
            continue
        gene_id = gene_sum.strip('.tsv.gz')
        df = pd.read_table(os.path.join(srd_dir, gene_sum))
        df['var_id'] = gene_id + '_chr' + df['phe_chrom'].astype(str) + '_' + df['var_from'].astype(str) + '_' + df[
            'ref'].astype(str) + '_' + df['alt'].astype(str) + '_b38'
        df_list.append(df)
    result_df = pd.concat(df_list)
    result_df.to_csv(dest_file, sep='\t', header=True, index=False)


if __name__ == '__main__':
    # convert_sbams_dosage('/Users/haiyue.meng/Downloads/sim_scheme/sim_data/gwas_sbams',
    #                      '/Users/haiyue.meng/Downloads/sem/gwas',
    #                      '/Users/haiyue.meng/Downloads/sem/gtf_all_gene.tsv')
    # convert_sbams_dosage('/Users/haiyue.meng/Downloads/sim_scheme/sim_data/eqtl_sbams',
    #                      '/Users/haiyue.meng/Downloads/sem/eqtl',
    #                      '/Users/haiyue.meng/Downloads/sem/gtf_all_gene.tsv')

    # convert_sbams_dosage('/Users/haiyue.meng/Downloads/sim_scheme/sim_data/gwas_sbams',
    #                      '/Users/haiyue.meng/Downloads/sem/g',
    #                      '/Users/haiyue.meng/Downloads/sem/gtf_all_gene.tsv')
    # convert_sbams_dosage('/Users/haiyue.meng/Downloads/sim_scheme/sim_data/eqtl_sbams',
    #                      '/Users/haiyue.meng/Downloads/sem/e',
    #                      '/Users/haiyue.meng/Downloads/sem/gtf_all_gene.tsv')
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--src_sbam_dir', dest='src_sbam_dir',
    #                     help='Directory of source sbams, file name must be *.sbams.dat')
    # parser.add_argument('--dest_dir', dest='dest_dir',
    #                     help='Output directory')
    # parser.add_argument('--type', dest='type',
    #                     help='split or merge')
    # args = parser.parse_args()
    # print(f'Accepted args:\n {args}')
    # if args.type == 'split':
    #     convert_sbams_to_split(args.src_sbam_dir, args.dest_dir)
    # else:
    #     convert_sbams_to_merged(args.src_sbam_dir, args.dest_dir, 'eqtl')

    update_eqtl_annotation()
    merge_sum()
    merge_sum(srd_dir='/Volumes/HD/biodata/sem/eqtl/sum',
              dest_file='/Volumes/HD/biodata/sem/merged_eqtl_sum_stats.tsv.gz')
