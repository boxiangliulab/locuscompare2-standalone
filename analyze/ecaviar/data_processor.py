import json
import os.path
import sys
from pathlib import Path
import pandas as pd
import logging
import re

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import constants as const, coloc_utils as utils, global_data_process as gdp, config


class ECaviarDataProcessor:
    ECAVIAR_TOOL_NAME = 'ecaviar'
    output_zscore_column_name1 = 'merged_info'
    output_zscore_column_name2 = 'zscore'
    output_zscore_columns = [output_zscore_column_name1, output_zscore_column_name2]
    shell_command_plink_execute = 'plink -vcf {} --r2 --matrix --out {}'

    def __init__(self):
        logging.info('init ECaviarEQTLProcessor')

    def prepare(self, working_dir, gwas_cluster_dir, gwas_cluster_summary, eqtl_group_dir, eqtl_report, ref_vcf_dir,
                gwas_col_dict, eqtl_col_dict, population, gwas_sample_size, eqtl_sample_size, var_id_col_name):
        eqtl_candidate_df = pd.read_csv(eqtl_report, sep=const.column_spliter)
        gwas_candidate_df = pd.read_csv(gwas_cluster_summary, sep=const.column_spliter)
        gwas_cluster_snps_dict = self.__get_cluster_significant_snps_dict(gwas_candidate_df, gwas_col_dict)
        output_base_dir = f'{working_dir}/candidate'
        Path(output_base_dir).mkdir(exist_ok=True, parents=True)

        for gwas_cluster_file in os.listdir(gwas_cluster_dir):
            # cluster file name example: chr{chromosome}_{position}-chr{}.tsv
            result = re.match('^(.*).tsv$', gwas_cluster_file)
            if result and not gwas_cluster_file.startswith('.'):
                report_variant_id = result.group(1)
                # get chr{x}_{position}
                variant_id = report_variant_id.split('-')[0]
                # gwas文件中显著位点的positions
                gwas_significant_positions = gwas_cluster_snps_dict[variant_id]
                chromosome = utils.get_chromosome_number(report_variant_id.split('_')[0])
                gwas_cluster_df = pd.read_csv(f'{gwas_cluster_dir}/{gwas_cluster_file}', sep=const.column_spliter)
                for eqtl_file in os.listdir(f'{eqtl_group_dir}/{chromosome}'):
                    gene_id = eqtl_file.split('.')[0]
                    # eqtl一个基因文件中显著位点的positions
                    eqtl_significant_positions = self.__get_eqtl_significant_snp_positions(eqtl_candidate_df, eqtl_file)
                    # 如果当前基因文件中至少有一个显著位点存在于当前gwas的显著位点里，才进行gwas和eqtl的位点mapping工作
                    if any(pos in gwas_significant_positions for pos in eqtl_significant_positions):
                        if gene_id.upper().startswith(const.gene_id_prefix):
                            eqtl_grouped_df = pd.read_csv(f'{eqtl_group_dir}/{chromosome}/{eqtl_file}',
                                                          sep=const.column_spliter)
                            # gwas, eqtl position matching的交集
                            matching_gwas_df = gwas_cluster_df[
                                gwas_cluster_df[gwas_col_dict['position']].isin(
                                    eqtl_grouped_df[eqtl_col_dict['position']])]

                            if len(matching_gwas_df) == 0:
                                continue

                            vcf_output_dir = f'{working_dir}/vcf'
                            Path(vcf_output_dir).mkdir(parents=True, exist_ok=True)
                            candidate_dir = f'{output_base_dir}/{gene_id.upper()}/{variant_id}'
                            Path(candidate_dir).mkdir(parents=True, exist_ok=True)
                            # 1.通过gwas, eqtl, vcf交集数据生成vcf, run plink 生成LD file
                            input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chromosome}.vcf.gz')
                            output_vcf_name = f'{variant_id}.vcf'
                            if Path(input_vcf).exists():
                                snp_varids_positions = matching_gwas_df[[var_id_col_name,
                                                                        gwas_col_dict['position']]]
                                snp_varids_positions.sort_values(by=gwas_col_dict['position'], inplace=True)
                                utils.extract_vcf_data(chromosome, snp_varids_positions, input_vcf, vcf_output_dir,
                                                       output_vcf_name,
                                                       gwas_col_dict['position'],
                                                       target_snp_col_name=var_id_col_name)

                                vcf_matching_file = os.path.join(vcf_output_dir, 'matching', f'{variant_id}.tsv')
                                if not os.path.exists(vcf_matching_file) or os.path.getsize(vcf_matching_file) <= 0:
                                    logging.warning(f'No generated vcf file for gene {gene_id}')
                                    continue
                                vcf_matching_df = pd.read_table(vcf_matching_file,
                                                                sep=const.column_spliter, header=0,
                                                                usecols=[gwas_col_dict['position'], var_id_col_name],
                                                                dtype={gwas_col_dict['position']: 'Int64'})

                                # vcf_matching_file 包含gwas, eqtl, vcf的交集snps
                                vcf_matching_df.sort_values(by=gwas_col_dict['position'], inplace=True)
                                output_ld_file = f'{candidate_dir}/{variant_id}_{gene_id}'
                                os.system(self.shell_command_plink_execute.format(f'{vcf_output_dir}/{output_vcf_name}',
                                                                                  output_ld_file))
                                # 移除LD值为nan的行和列，并返回被移除的snp positions
                                nan_cols = utils.remove_nan_from_ld(f'{output_ld_file}.ld',
                                                                    vcf_matching_df[gwas_col_dict['position']].tolist())

                                # 从vcf matching snps中移除LD行列为nan的snp
                                vcf_matching_df = vcf_matching_df[~vcf_matching_df[gwas_col_dict['position']].isin(nan_cols)]

                                # Drop GWAS rows that does not have vcf records
                                # SNP in vcf_matching_df is subset of SNP in candidate_gwas_df, so it's fine to drop intersect rows here
                                utils.drop_non_intersect_rows(matching_gwas_df, var_id_col_name, vcf_matching_df,
                                                              var_id_col_name)
                                # Drop eQTL rows that does not have vcf records
                                utils.drop_non_intersect_rows(eqtl_grouped_df, var_id_col_name, vcf_matching_df,
                                                              var_id_col_name)

                                # clean vcf file after ld computation
                                # utils.delete_file_if_exists(f'{vcf_output_dir}/{output_vcf_name}')

                            else:
                                logging.warning(f'!ref vcf {input_vcf} does not exist')

                            # 2.根据gwas, eqtl, vcf position交集数据生成eqtl zscore
                            eqtl_grouped_df.sort_values(by=eqtl_col_dict['position'], inplace=True)
                            eqtl_grouped_df[self.output_zscore_column_name1] = eqtl_grouped_df[
                                                                                    eqtl_col_dict['chrom']].map('{}:'.format) \
                                                                                + eqtl_grouped_df[
                                                                                    eqtl_col_dict['position']].apply(str) \
                                                                                + eqtl_grouped_df[
                                                                                    eqtl_col_dict['alt']].map(':{}:'.format) \
                                                                                + eqtl_grouped_df[eqtl_col_dict['ref']]
                            eqtl_grouped_df[self.output_zscore_column_name2] = eqtl_grouped_df[
                                                                                    eqtl_col_dict['beta']] / \
                                                                                eqtl_grouped_df[eqtl_col_dict['se']]
                            eqtl_zscore_file = f'{candidate_dir}/eqtl_{variant_id}_{gene_id}.z'
                            eqtl_grouped_df.to_csv(eqtl_zscore_file, columns=self.output_zscore_columns, header=None,
                                                    sep=' ', index=False)

                            # 3.根据gwas, eqtl, vcf position交集数据生成gwas zscore
                            utils.adjust_allele_order(matching_gwas_df,
                                                      gwas_col_dict['effect_allele'],
                                                      gwas_col_dict['other_allele'],
                                                      gwas_col_dict['chrom'],
                                                      gwas_col_dict['position'],
                                                      eqtl_grouped_df,
                                                      ref_df_chrom_col_name=eqtl_col_dict['chrom'],
                                                      ref_df_pos_col_name=eqtl_col_dict['position'],
                                                      ref_df_alt_allele_col_name=eqtl_col_dict['alt'],
                                                      ref_df_ref_allele_col_name=eqtl_col_dict['ref'],
                                                      gbeta_col_name=gwas_col_dict['beta'])
                            matching_gwas_df.sort_values(by=gwas_col_dict['position'], inplace=True)
                            matching_gwas_df[self.output_zscore_column_name1] = matching_gwas_df[
                                                                                    gwas_col_dict['chrom']].map('{}:'.format) \
                                                                                + matching_gwas_df[
                                                                                    gwas_col_dict['position']].apply(str) \
                                                                                + matching_gwas_df[
                                                                                    gwas_col_dict['effect_allele']].map(':{}:'.format) \
                                                                                + matching_gwas_df[
                                                                                    gwas_col_dict['other_allele']]
                            matching_gwas_df[self.output_zscore_column_name2] = matching_gwas_df[
                                                                                    gwas_col_dict['beta']] / \
                                                                                matching_gwas_df[gwas_col_dict['se']]
                            gwas_zscore_file = f'{candidate_dir}/gwas_{variant_id}_{gene_id}.z'
                            matching_gwas_df.to_csv(gwas_zscore_file, columns=self.output_zscore_columns, header=None,
                                                    sep=' ', index=False)

                            # 生成finemap所需的in file
                            finemap_in_file = f'{candidate_dir}/{variant_id}_{gene_id}.in'
                            with open(finemap_in_file, mode='w') as in_file:
                                in_file.write('z;ld;snp;config;n-ind\n')
                                in_file.write(f'{gwas_zscore_file};{output_ld_file}.ld;{candidate_dir}/gwas_{variant_id}_{gene_id}.snp;{candidate_dir}/gwas_{variant_id}_{gene_id}.config;{gwas_sample_size}\n')
                                in_file.write(f'{eqtl_zscore_file};{output_ld_file}.ld;{candidate_dir}/eqtl_{variant_id}_{gene_id}.snp;{candidate_dir}/eqtl_{variant_id}_{gene_id}.config;{eqtl_sample_size}')
        return output_base_dir

    def __get_cluster_significant_snps_dict(self, cluster_df, gwas_col):
        grouped = cluster_df.groupby('range_lead')
        cluster_snps_dict = {}
        for name, group in grouped:
            cluster_snps_dict[name] = group[gwas_col['position']].tolist()
        return cluster_snps_dict

    def __get_eqtl_significant_snp_positions(self, eqtl_df, gene_file_name):
        results = eqtl_df[eqtl_df['gene_file'] == gene_file_name]
        if results.empty:
            return []
        else:
            significant_positions = results.iloc[0]['positions']
            if isinstance(significant_positions, str):
                return json.loads(significant_positions)
            elif isinstance(significant_positions, list):
                return significant_positions
            else:
                return []


if __name__ == '__main__':
    cfg = config.ConfigHolder()
    glob_processor = gdp.Processor(cfg_holder=cfg)
    _working_dir = os.path.join('/Users/nicklin/Desktop/CAD/test', ECaviarDataProcessor.ECAVIAR_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    pop = 'EUR'
    gwas_data_processor = ECaviarDataProcessor()

    _gwas_col_dict = {'snp': 'oldID', 'chrom': 'chr', 'position': 'hm_pos', 'beta': 'beta',
                      'effect_allele': 'a1', 'other_allele': 'a2', 'pvalue': 'pval',
                      'se': 'se', 'eaf': 'eaf'}
    _eqtl_col_dict = {'snp': 'rsid', 'chrom': 'chromosome', 'position': 'position', 'beta': 'beta', 'alt': 'alt',
                      'ref': 'ref', 'pvalue': 'pvalue', 'se': 'se', 'gene_id': 'molecular_trait_id', 'maf': 'maf'}

    gwas_data_processor.prepare(working_dir=_working_dir,
                                gwas_cluster_dir='/Volumes/HD/biodata/colocalization-tools/preprocessed/gwas/cad/clustered',
                                gwas_cluster_summary='/Volumes/HD/biodata/colocalization-tools/preprocessed/gwas/cad/cluster_summary.tsv',
                                eqtl_group_dir='/Volumes/HD/biodata/colocalization-tools/preprocessed/eqtl/Artery_Aorta/grouped',
                                eqtl_report='/Volumes/HD/biodata/colocalization-tools/preprocessed/eqtl/Artery_Aorta/filtered_gene.tsv',
                                ref_vcf_dir='/Volumes/HD/biodata/colocalization-tools/raw/vcf/hg38_phased',
                                gwas_col_dict=_gwas_col_dict,
                                eqtl_col_dict=_eqtl_col_dict,
                                population=pop,
                                gwas_sample_size=296525,
                                eqtl_sample_size=948,
                                var_id_col_name='var_id')
