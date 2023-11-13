import concurrent
import datetime
import logging
import os
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.append(os.path.abspath(os.path.dirname(Path(__file__).resolve())))
import coloc_utils as utils
import constants as const

'''
预处理的GWAS文件路径:
{working_dir}/preprocessed/gwas/{trait}
预处理的eQTL文件路径:
{working_dir}/preprocessed/eqtl/{tissue}

工具working dir:
{working_dir}/processed/{trait}_{tissue}/{tool}
{working_dir}/processed/{trait}_{tissue}/{tool}
'''


class Processor:
    # column to identify a unique record in gwas & eqtl files, in chromosome_position format
    VAR_ID_COL_NAME = 'var_id_'

    def __init__(self, cfg_holder=None):
        logging.info('init Processor')
        self.config_holder = cfg_holder
        self.global_config = self.config_holder.global_config
        self.root_work_dir = self.config_holder.root_work_dir
        self.output_preprocessed_dir = self.config_holder.output_preprocessed_dir
        self.output_processed_dir = self.config_holder.output_processed_dir

        self.gwas_trait = self.config_holder.gwas_trait
        self.eqtl_tissue = self.config_holder.eqtl_tissue

        self.tool_parent_dir = self.config_holder.tool_parent_dir
        self.rank_dir = self.config_holder.rank_dir
        self.report_path = self.config_holder.report_path

        self.gwas_col_dict = self.config_holder.gwas_col_dict
        self.eqtl_col_dict = self.config_holder.eqtl_col_dict

        # eqtl output path
        self.eqtl_output_dir = self.config_holder.eqtl_output_dir
        self.eqtl_output_report = self.config_holder.eqtl_output_report

        # gwas output path
        self.gwas_preprocessed_dir = self.config_holder.gwas_preprocessed_dir
        self.gwas_preprocessed_file = self.config_holder.gwas_preprocessed_file
        self.gwas_filter_file = self.config_holder.gwas_filter_file
        self.gwas_cluster_output_dir = self.config_holder.gwas_cluster_output_dir
        self.gwas_output_dir = self.config_holder.gwas_output_dir
        self.gwas_cluster_summary = self.config_holder.gwas_cluster_summary

        # vcf input output path
        self.vcf_output_dir = self.config_holder.vcf_output_dir
        self.ref_vcf_dir = self.config_holder.ref_vcf_dir

        # tools parameter config path
        self.tools_config_file = self.config_holder.tools_config_file

    def preprocess_eqtl(self):
        # group by trait
        logging.info('start to split eQTL file by gene')
        utils.delete_dir(self.eqtl_output_dir)
        Path(self.eqtl_output_dir).mkdir(exist_ok=True, parents=True)
        utils.split_file_by_col_name(self.eqtl_output_dir, f'{self.global_config["input"]["eqtl"]["file"]}',
                                     self.eqtl_col_dict['chrom'], self.eqtl_col_dict['gene_id'],
                                     readonly_cols=self.eqtl_col_dict.values(), sep=self.config_holder.eqtl_sep)
        logging.info('finish split eQTL file')
        # drop na or no significant snp gene file, sort by position
        total_gene_file_count = 0
        good_result_count = 0
        chrom_list = []
        gene_file_list = []
        postions_list = []
        logging.info('start to filter gene eQTL file')
        logging.info(f'Filtering eQTL data by p-value threshold {self.config_holder.eqtl_p_threshold}')
        for chrom_dir in os.listdir(f'{self.eqtl_output_dir}'):
            if not chrom_dir.startswith('.'):
                for eqtl_file in os.listdir(f'{self.eqtl_output_dir}/{chrom_dir}'):
                    if eqtl_file.upper().startswith(const.gene_id_prefix):
                        total_gene_file_count += 1
                        _, high_risk_snp_df = self.__prepare_eqtl_data(
                            f'{self.eqtl_output_dir}/{chrom_dir}/{eqtl_file}')
                        if high_risk_snp_df is not None:
                            good_result_count += 1
                            chrom_list.append(chrom_dir)
                            gene_file_list.append(eqtl_file)
                            postions_list.append(high_risk_snp_df[self.eqtl_col_dict['position']].tolist())
                            logging.debug(f'good gene {eqtl_file}')
        eqtl_filter_result = pd.DataFrame(
            {'chrom': chrom_list, 'gene_file': gene_file_list, 'positions': postions_list})
        eqtl_filter_result.to_csv(self.eqtl_output_report, sep=const.output_spliter, index=False)
        logging.info(
            f'finish filter gene eQTL file, {total_gene_file_count} gene files, good result {good_result_count}')
        return eqtl_filter_result

    def __prepare_eqtl_data(self, eqtl_trait_file_path):
        try:
            eqtl_trait_df = pd.read_table(eqtl_trait_file_path, sep=const.column_spliter, header=0,
                                          usecols=self.eqtl_col_dict.values(),
                                          dtype={self.eqtl_col_dict['chrom']: 'category',
                                                 self.eqtl_col_dict['position']: 'Int64'})
        except Exception as e:
            logging.error(f'error to prepare eqtl data: {eqtl_trait_file_path}')
            logging.error(f'exception: {e}')
        eqtl_trait_df[Processor.VAR_ID_COL_NAME] = 'chr' + eqtl_trait_df[self.eqtl_col_dict['chrom']].astype(
            str) + '_' + eqtl_trait_df[self.eqtl_col_dict['position']].astype(str)
        utils.clean_data(eqtl_trait_df, dup_consider_subset=Processor.VAR_ID_COL_NAME)
        utils.drop_indel_snp(eqtl_trait_df, self.eqtl_col_dict['alt'], self.eqtl_col_dict['ref'])
        if eqtl_trait_df.empty:
            del eqtl_trait_df
            utils.delete_file_if_exists(eqtl_trait_file_path)
            logging.info(f'No required data in eQTL file {eqtl_trait_file_path}')
            return None, None
        else:
            pval_filter_eqtl_df = \
                utils.filter_data_frame_by_p_value(eqtl_trait_df, self.config_holder.eqtl_p_threshold,
                                                   self.eqtl_col_dict['pvalue'], inplace=False)
            if pval_filter_eqtl_df.empty:
                del eqtl_trait_df
                del pval_filter_eqtl_df
                utils.delete_file_if_exists(eqtl_trait_file_path)
                return None, None
            else:
                eqtl_trait_df.sort_values([self.eqtl_col_dict['chrom'], self.eqtl_col_dict['position']],
                                          inplace=True)
                eqtl_trait_df.to_csv(eqtl_trait_file_path, sep=const.output_spliter, index=False)
                return eqtl_trait_df, pval_filter_eqtl_df

    def __merge_alt_ref(self, gwas_df, chrom, population):
        logging.info(f'Merging alt ref for chrom {chrom}')
        input_vcf = os.path.join(self.ref_vcf_dir, population, f'chr{chrom}.vcf.gz')
        vcf_df = pd.read_table(input_vcf, header=None, comment='#', usecols=[0, 1, 3, 4],
                               dtype={0: 'category', 1: 'Int64', 3: pd.CategoricalDtype(const.SNP_ALLELE),
                                      4: pd.CategoricalDtype(const.SNP_ALLELE)})
        vcf_df.columns = [self.gwas_col_dict['chrom'], self.gwas_col_dict['position'], 'ref_', 'alt_']
        vcf_df.dropna(subset=['ref_', 'alt_'], inplace=True)
        vcf_df.reset_index(drop=True, inplace=True)
        merged_vcf_col_suffix = f'_{chrom}_vcf'
        gwas_df = pd.merge(left=gwas_df,
                           right=vcf_df,
                           left_on=[self.gwas_col_dict['chrom'], self.gwas_col_dict['position']],
                           right_on=[self.gwas_col_dict['chrom'], self.gwas_col_dict['position']],
                           how='left',
                           suffixes=(None, merged_vcf_col_suffix))
        del vcf_df
        gwas_df['ref_'].mask((gwas_df[self.gwas_col_dict['chrom']] == chrom) & gwas_df['ref_'].isna(),
                             gwas_df[f'ref_{merged_vcf_col_suffix}'], inplace=True)
        gwas_df['alt_'].mask((gwas_df[self.gwas_col_dict['chrom']] == chrom) & gwas_df['alt_'].isna(),
                             gwas_df[f'alt_{merged_vcf_col_suffix}'], inplace=True)
        gwas_df.drop(columns=[col for col in gwas_df.columns if col.endswith(merged_vcf_col_suffix)], inplace=True)
        logging.info(f'Merge alt ref for chrom {chrom} completed')

    # def __gwas_clumping(self, gwas_chrom_file, chrom, population):
    #     logging.info(f'Clumping for chrom {chrom}')
    #     input_vcf = os.path.join(self.ref_vcf_dir, population, f'chr{chrom}.vcf.gz')
    #     clumped_out = os.path.join(self.gwas_output_dir, f'chr{chrom}')
    #     clump_r2 = self.global_config.get('clump_r2', 0)
    #     clump_kb = self.global_config.get('clump_kb', 500)
    #     logging.info(f'Clumping param for chrom {chrom} clump_r2: {clump_r2}, clump_kb: {clump_kb}')
    #     os.system(f'plink --silent --vcf {input_vcf}  '
    #               f'--clump {gwas_chrom_file} '
    #               f'--clump-p1 {self.config_holder.gwas_p_threshold} '
    #               f'--clump-p2 {self.config_holder.gwas_p_threshold} '
    #               f'--clump-snp-field {self.gwas_col_dict["snp"]} '
    #               f'--clump-field {self.gwas_col_dict["pvalue"]} '
    #               f'--clump-r2 {clump_r2} '
    #               f'--clump-kb {clump_kb} '
    #               f'--out {clumped_out}')
    #     logging.info(f'Clumping for chrom {chrom} completed')

    def preprocess_gwas(self):
        utils.delete_dir(self.gwas_preprocessed_dir)
        Path(self.gwas_preprocessed_dir).mkdir(exist_ok=True, parents=True)
        gwas_file_path = self.global_config['input']['gwas']['file']
        population = self.global_config.get('population', 'EUR').upper()
        logging.info(f'Reading GWAS file {gwas_file_path}')
        gwas_df = pd.read_table(gwas_file_path, sep=self.config_holder.gwas_sep, header=0,
                                usecols=self.gwas_col_dict.values(),
                                dtype={self.gwas_col_dict['position']: 'Int64', self.gwas_col_dict['chrom']: 'category',
                                       self.gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                       self.gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)})
        raw_gwas_snp_size = len(gwas_df)
        logging.info(f'GWAS data memory usage: {gwas_df.memory_usage(deep=True)}')
        gwas_df[self.gwas_col_dict['chrom']] = gwas_df[self.gwas_col_dict['chrom']].astype(str).str.lower().str.strip(
            'chr').astype('category')
        logging.info(f'GWAS data dropping non-autosome data, time: {datetime.datetime.now()}')
        gwas_df.drop(labels=gwas_df[~gwas_df[self.gwas_col_dict['chrom']].isin([str(i) for i in range(1, 23)])].index,
                     inplace=True)
        logging.info(f'GWAS data dropping INDEL SNPs, time: {datetime.datetime.now()}')
        utils.drop_indel_snp(gwas_df, self.gwas_col_dict['effect_allele'], self.gwas_col_dict['other_allele'])
        logging.info(f'GWAS data cleaning, time: {datetime.datetime.now()}')
        utils.clean_data(gwas_df, dup_consider_subset=[self.gwas_col_dict['chrom'], self.gwas_col_dict['position']])
        gwas_df.drop(index=gwas_df[gwas_df[self.gwas_col_dict['se']] == 0].index, inplace=True)
        filtered_gwas_snp_size = len(gwas_df)
        logging.info(
            f'GWAS data total {filtered_gwas_snp_size} rows, sorting by chrom and pos, time: {datetime.datetime.now()}')
        gwas_df.sort_values([self.gwas_col_dict['chrom'], self.gwas_col_dict['position']], inplace=True)
        logging.info(f'Merging alt ref from vcf into GWAS data')
        # discontinuous index cost a lot more memory
        gwas_df.reset_index(drop=True, inplace=True)
        # Merge alt/ref from vcf into gwas file for later use
        gwas_df['alt_'] = pd.Series(data=pd.NA, dtype=pd.CategoricalDtype(const.SNP_ALLELE))
        gwas_df['ref_'] = pd.Series(data=pd.NA, dtype=pd.CategoricalDtype(const.SNP_ALLELE))
        chroms = gwas_df[self.gwas_col_dict['chrom']].unique()
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = []
            for chrom in chroms:
                futures.append(executor.submit(self.__merge_alt_ref, gwas_df, chrom, population))
            for future in concurrent.futures.as_completed(futures):
                try:
                    data = future.result()
                except Exception as exc:
                    logging.error("".join(traceback.TracebackException.from_exception(exc).format()))
        gwas_df.drop_duplicates(subset=[self.gwas_col_dict['chrom'], self.gwas_col_dict['position']],
                                keep=False, inplace=True)
        gwas_df[self.gwas_col_dict['effect_allele']] = gwas_df[self.gwas_col_dict['effect_allele']].str.upper().astype(
            pd.CategoricalDtype(const.SNP_ALLELE))
        gwas_df[self.gwas_col_dict['other_allele']] = gwas_df[self.gwas_col_dict['other_allele']].str.upper().astype(
            pd.CategoricalDtype(const.SNP_ALLELE))
        # Fill SNP ref/alt by other_allele/effect_allele if the SNP is not in vcf
        gwas_df['ref_'].mask(gwas_df['ref_'].isna(), gwas_df[self.gwas_col_dict['other_allele']], inplace=True)
        gwas_df['alt_'].mask(gwas_df['alt_'].isna(), gwas_df[self.gwas_col_dict['effect_allele']], inplace=True)
        gwas_df[Processor.VAR_ID_COL_NAME] = 'chr' + gwas_df[self.gwas_col_dict['chrom']].astype(str) + '_' + \
                                             gwas_df[self.gwas_col_dict['position']].astype(str)
        logging.info(
            f'Writing GWAS preprocessed data to {self.gwas_preprocessed_file}, time: {datetime.datetime.now()}')
        gwas_df.to_csv(self.gwas_preprocessed_file, sep=const.output_spliter, header=True, index=False)
        logging.info(f'Filtering GWAS data by p-value threshold {self.config_holder.gwas_p_threshold}')
        pval_filter_gwas_df = \
            utils.filter_data_frame_by_p_value(gwas_df,
                                               self.config_holder.gwas_p_threshold,
                                               self.gwas_col_dict['pvalue'],
                                               inplace=False)
        logging.info(
            f'GWAS data filtered result {len(pval_filter_gwas_df)} rows, writing to file, time: {datetime.datetime.now()}')
        pval_filter_gwas_df.to_csv(self.gwas_filter_file, sep=const.output_spliter, header=True, index=False)
        if pval_filter_gwas_df.empty:
            # print(f'No significant records found in GWAS file {gwas_file_path}')
            logging.warning(f'No significant records found in GWAS file {gwas_file_path}')
            return
        utils.delete_dir(self.gwas_output_dir)
        Path(self.gwas_output_dir).mkdir(exist_ok=True, parents=True)
        utils.delete_dir(self.gwas_cluster_output_dir)
        Path(self.gwas_cluster_output_dir).mkdir(exist_ok=True, parents=True)
        pval_filter_gwas_df.sort_values([self.gwas_col_dict['chrom'], self.gwas_col_dict['pvalue']], inplace=True)
        chrom_list = []
        range_lead_list = []
        positions_list = []
        cluster_pos_dict = {}
        loci_count = 0
        clump_dist = 500000
        for _, row in pval_filter_gwas_df.iterrows():
            chrom = row.loc[self.gwas_col_dict['chrom']]
            pos = row.loc[self.gwas_col_dict['position']]
            chrom_cluster_pos_list = cluster_pos_dict.get(chrom, [])
            for cp in chrom_cluster_pos_list:
                cp_start = max(0, cp - clump_dist)
                cp_end = cp + clump_dist
                if cp_start <= pos <= cp_end:
                    break
            else:
                chrom_cluster_pos_list.append(pos)
                cluster_pos_dict[chrom] = chrom_cluster_pos_list
                # retrieve range_df from group by peak_positions.min <= group[position] <= peak_positions.max
                cluster_start = max(0, pos - clump_dist)
                cluster_end = pos + clump_dist
                range_df = gwas_df[(gwas_df[self.gwas_col_dict['chrom']] == chrom) & (
                        cluster_start <= gwas_df[self.gwas_col_dict['position']]) & (
                                           gwas_df[self.gwas_col_dict['position']] <= cluster_end)]
                range_lead_var_id = \
                    range_df.loc[range_df[self.gwas_col_dict['position']] == pos][Processor.VAR_ID_COL_NAME].iloc[0]
                chrom_list.append(chrom)
                range_lead_list.append(range_lead_var_id)
                loci_count += 1
                file_name = f'{range_lead_var_id}-chr{chrom}.tsv.gz'
                range_df.to_csv(os.path.join(self.gwas_cluster_output_dir, file_name), sep=const.output_spliter,
                                header=True,
                                index=False)
                del range_df
                pval_filter_range_df = pval_filter_gwas_df[(pval_filter_gwas_df[self.gwas_col_dict['chrom']] == chrom) & (
                        cluster_start <= pval_filter_gwas_df[self.gwas_col_dict['position']]) & (
                                                                   pval_filter_gwas_df[
                                                                       self.gwas_col_dict['position']] <= cluster_end)]
                positions_list.append(pval_filter_range_df[self.gwas_col_dict['position']].tolist())
                del pval_filter_range_df
        del pval_filter_gwas_df
        # neighbor_range = self.global_config['neighbour_snp_range']
        logging.info(
            f'Clumping GWAS SNPs, range files will be written to {self.gwas_cluster_output_dir}, time: {datetime.datetime.now()}')
        gwas_chrom_group_files = {}
        for name, group in gwas_df.groupby(self.gwas_col_dict['chrom'], observed=True):
            if group.shape[0] == 0:
                continue
            group_file = os.path.join(self.gwas_output_dir, f'chr{name}.tsv.gz')
            group.to_csv(group_file, sep=const.output_spliter, header=True, index=False)
            gwas_chrom_group_files[name] = group_file
        del gwas_df
        # with ThreadPoolExecutor(max_workers=4) as executor:
        #     futures = []
        #     for name, group_file in gwas_chrom_group_files.items():
        #         futures.append(executor.submit(self.__gwas_clumping, group_file, name, population))
        #     for future in concurrent.futures.as_completed(futures):
        #         try:
        #             data = future.result()
        #         except Exception as exc:
        #             logging.error('Get result generated an exception: %s' % exc)
        # loci_count = 0
        # for name, group_file in gwas_chrom_group_files.items():
        #     group = pd.read_table(group_file, sep=const.output_spliter,
        #                           dtype={self.gwas_col_dict['position']: 'Int64',
        #                                  self.gwas_col_dict['chrom']: 'category',
        #                                  self.gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
        #                                  self.gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)})
        #     clumped_out = os.path.join(self.gwas_output_dir, f'chr{name}')
        #     clumped_file = f'{clumped_out}.clumped'
        #     if not os.path.exists(clumped_file) or os.path.getsize(clumped_file) <= 0:
        #         logging.warning(f'No significant SNP on chromosome {name}')
        #         continue
        #     index_variant_df = pd.read_table(clumped_file, sep=r'\s+', usecols=['SNP', 'BP', 'TOTAL', 'SP2'])
        #     index_variant_df.drop(index=index_variant_df[index_variant_df['TOTAL'] == 0].index, inplace=True)
        #     loci_count += index_variant_df.shape[0]
        #     for _, row in index_variant_df.iterrows():
        #         peak_snps = re.sub(r'\(\d+\)', '', row.loc['SP2']).split(',')
        #         range_df = utils.find_peak_range_df(name, row.loc['BP'], group, self.gwas_col_dict['chrom'],
        #                                             self.gwas_col_dict['position'], neighbor_range)
        #         if len(peak_snps) == 0 or peak_snps[0] == 'NONE':
        #             peak_positions = [row.loc['BP']]
        #         else:
        #             # consider the index snp(i.e. lead snp), index snp maybe on the edge of the cluster
        #             peak_snps.append(row.loc['SNP'])
        #             # NOTE!! The column gwas_col_dict['snp'] in gwas file can be other values(like variant_id),
        #             # as long as they match the values of ID column in vcf file, else clumping won't work
        #             peak_positions = group.loc[group[self.gwas_col_dict['snp']].isin(peak_snps)][
        #                 self.gwas_col_dict['position']].tolist()
        #         # all peak_positions should be in range_df, else extends range_df
        #         if not all(pos in range_df[self.gwas_col_dict['position']].values for pos in peak_positions):
        #             # retrieve range_df from group by peak_positions.min <= group[position] <= peak_positions.max
        #             range_df = group[(min(peak_positions) <= group[self.gwas_col_dict['position']]) & (
        #                     group[self.gwas_col_dict['position']] <= max(peak_positions))]
        #         range_lead_var_id = \
        #             range_df.loc[range_df[self.gwas_col_dict['position']] == row.loc['BP']][
        #                 Processor.VAR_ID_COL_NAME].iloc[0]
        #         positions_list.append(peak_positions)
        #         chrom_list.append(name)
        #         range_lead_list.append(range_lead_var_id)
        #         file_name = f'{range_lead_var_id}-chr{name}.tsv.gz'
        #         range_df.to_csv(os.path.join(self.gwas_cluster_output_dir, file_name), sep=const.output_spliter,
        #                         header=True,
        #                         index=False)
        #     del group
        with open(f'{self.gwas_preprocessed_dir}/gwas_snp_summary.log', mode='w') as snp_file:
            snp_file.write(f'raw GWAS snp size: {raw_gwas_snp_size}\n')
            snp_file.write(f'filtered GWAS snp size: {filtered_gwas_snp_size}\n')
            snp_file.write(f'GWAS significant loci size: {loci_count}\n')
        logging.info(
            f'Writing GWAS significant ranges summary to {self.gwas_cluster_summary}, time: {datetime.datetime.now()}')
        cluster_summary_df = pd.DataFrame(
            {'chrom': chrom_list, 'range_lead': range_lead_list, 'positions': positions_list})
        cluster_summary_df.to_csv(self.gwas_cluster_summary, sep=const.output_spliter, header=True, index=False)
        logging.info(f'GWAS data preprocess completed, time: {datetime.datetime.now()}')


if __name__ == '__main__':
    processor = Processor()
    processor.preprocess_gwas()
    processor.preprocess_eqtl()
