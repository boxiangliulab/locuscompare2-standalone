import concurrent
import datetime
import logging
import os
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import ast

import pandas as pd

sys.path.append(os.path.abspath(os.path.dirname(Path(__file__).resolve())))
import utils
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
        self.qtl_type = self.config_holder.qtl_type
        self.gwas_trait = self.config_holder.gwas_trait
        self.biological_context = self.config_holder.biological_context

        self.gwas_genomic_window = self.config_holder.gwas_genomic_window
        self.gwas_window_size = self.config_holder.gwas_window_size
        self.gwas_LD_r2_filter = self.config_holder.gwas_LD_r2_filter
        self.gwas_LD_additional_expansion = self.config_holder.gwas_LD_additional_expansion

        self.qtl_genomic_window = self.config_holder.qtl_genomic_window
        self.qtl_window_size = self.config_holder.qtl_window_size
        self.qtl_LD_r2_filter = self.config_holder.qtl_LD_r2_filter
        self.qtl_LD_additional_expansion = self.config_holder.qtl_LD_additional_expansion

        self.global_LD_based_window = self.config_holder.global_LD_based_window
        print(f"self.global_LD_based_window : {self.global_LD_based_window}")

        self.target_loci = self.config_holder.target_loci

        self.tool_parent_dir = self.config_holder.tool_parent_dir
        self.rank_dir = self.config_holder.rank_dir
        self.report_path = self.config_holder.report_path

        self.gwas_col_dict = self.config_holder.gwas_col_dict
        self.qtl_col_dict = self.config_holder.qtl_col_dict

        # qtl output path
        self.qtl_preprocesed_dir = self.config_holder.qtl_preprocesed_dir
        self.qtl_grouped_dir = self.config_holder.qtl_grouped_dir
        self.qtl_output_report = self.config_holder.qtl_output_report
        self.qtl_LD_window = self.config_holder.qtl_LD_window

        # gwas output path
        self.gwas_preprocessed_dir = self.config_holder.gwas_preprocessed_dir
        self.gwas_preprocessed_file = self.config_holder.gwas_preprocessed_file
        self.gwas_filter_file = self.config_holder.gwas_filter_file
        self.gwas_cluster_output_dir = self.config_holder.gwas_cluster_output_dir
        self.gwas_output_dir = self.config_holder.gwas_output_dir
        self.gwas_cluster_summary = self.config_holder.gwas_cluster_summary # preprocessed/gwas/trait/cluster_summary.tsv.gz

        # vcf input output path
        self.vcf_output_dir = self.config_holder.vcf_output_dir
        self.ref_vcf_dir = self.config_holder.ref_vcf_dir
        self.min_matching_number = self.config_holder.min_matching_number

        # tools parameter config path
        self.tools_config_file = self.config_holder.tools_config_file

        self.shell_command_plink_execute = 'plink --silent --vcf {} --r2 --ld-snp {} --ld-window-kb {} --ld-window-r2 {} --out {}'


    def __merge_alt_ref(self, gwas_df, chrom, population):        
        logging.info(f'Merging alt ref for chrom {chrom}')
        input_vcf = os.path.join(self.ref_vcf_dir, population, f'chr{chrom}.vcf.gz')
        vcf_df = pd.read_table(input_vcf, header=None, comment='#', usecols=[0, 1, 3, 4],
                               dtype={0: 'category', 
                                      1: 'Int64', 
                                    #   3: pd.CategoricalDtype(const.SNP_ALLELE),
                                    #   4: pd.CategoricalDtype(const.SNP_ALLELE)
                                      3: 'category',
                                      4: 'category'
                                      })
        vcf_df.columns = [self.gwas_col_dict['chrom'], self.gwas_col_dict['position'], 'ref_', 'alt_']
        vcf_df.dropna(subset=['ref_', 'alt_'], inplace=True)
        vcf_df.reset_index(drop=True, inplace=True)
        merged_vcf_col_suffix = f'_{chrom}_vcf'
        print(f"merged_vcf_col_suffix: {merged_vcf_col_suffix}")
        # print(f"__merge_alt_ref gwas_df: {gwas_df}")
        vcf_df['ref_'] = vcf_df['ref_'].str.upper()
        vcf_df['alt_'] = vcf_df['alt_'].str.upper()
        vcf_df.index = vcf_df[self.gwas_col_dict['chrom']].astype(str) + '_' + vcf_df[self.gwas_col_dict['position']].astype(str) + '_' + vcf_df['ref_'] + '_' + vcf_df['alt_']
        gwas_df.index = gwas_df[self.gwas_col_dict['variant_id']]
        gwas_df = pd.merge(left=gwas_df,
                           right=vcf_df,
                           left_index=True,
                           right_index=True,
                           how='left',
                           suffixes=(None, merged_vcf_col_suffix))
        del vcf_df
        gwas_df['ref_'].mask((gwas_df[self.gwas_col_dict['chrom']] == chrom) & gwas_df['ref_'].isna(),
                             gwas_df[f'ref_{merged_vcf_col_suffix}'], inplace=True)
        gwas_df['alt_'].mask((gwas_df[self.gwas_col_dict['chrom']] == chrom) & gwas_df['alt_'].isna(),
                             gwas_df[f'alt_{merged_vcf_col_suffix}'], inplace=True)
        gwas_df.drop(columns=[col for col in gwas_df.columns if col.endswith(merged_vcf_col_suffix)], inplace=True)
        logging.info(f'Merge alt ref for chrom {chrom} completed')
################################################################################
#                             1.3 Preprocess eQTL                              #
################################################################################
#                    gene_id              variant_id  tss_distance        af  ma_samples  ma_count  pval_nominal     slope  slope_se
# 0        ENSG00000277117.5   chr21_5033246_G_T_b38         10753  0.016875          27        27  1.744312e-34 -1.889374  0.146365
# 1        ENSG00000277117.5   chr21_5033539_G_A_b38         11046  0.967500          52        52  7.627077e-54  1.649587  0.098137

    def preprocess_eqtl(self):
        # group by trait
        logging.info('Start to split eQTL file by gene')
        utils.delete_dir(self.qtl_grouped_dir)
        Path(self.qtl_grouped_dir).mkdir(exist_ok=True, parents=True)
        utils.split_file_by_col_name_qtl(self.qtl_grouped_dir, f'{self.global_config["input"]["qtl"]["file"]}',
                                     self.qtl_col_dict['chrom'], self.qtl_col_dict['phenotype_id'],
                                     readonly_cols=self.qtl_col_dict.values(), sep=self.config_holder.qtl_sep)
        logging.info('finish split eQTL file')
        # drop na or no significant snp gene file, sort by position
        total_gene_file_count = 0
        good_result_count = 0
        chrom_list = []
        gene_file_list = []
        postions_list = []
        logging.info('start to filter gene eQTL file')
        # logging.info(f'Filtering eQTL data by p-value threshold {self.config_holder.qtl_p_threshold}')
        for chrom_dir in os.listdir(f'{self.qtl_grouped_dir}'):
        # for chrom_dir in ['15','6']:
            if not chrom_dir.startswith('.'):
                for eqtl_file in os.listdir(f'{self.qtl_grouped_dir}/{chrom_dir}'):
                    # if eqtl_file.upper().startswith(const.gene_id_prefix):
                    total_gene_file_count += 1
                    gene_id = utils.get_eqtl_gene_name(eqtl_file)
                    high_risk_snp_df = self.__prepare_eqtl_data(
                        f'{self.qtl_grouped_dir}/{chrom_dir}/{eqtl_file}', gene_id)
                    if high_risk_snp_df is not None:
                        good_result_count += 1
                        chrom_list.append(chrom_dir)
                        gene_file_list.append(eqtl_file)
                        postions_list.append(high_risk_snp_df[self.qtl_col_dict['position']].tolist())
                        logging.debug(f'good gene {eqtl_file}')
        eqtl_filter_result = pd.DataFrame(
            {'chrom': chrom_list, 'pheno_file': gene_file_list, 'positions': postions_list})
        eqtl_filter_result.to_csv(self.qtl_output_report, sep=const.output_spliter, index=False)
        logging.info(
            f'finish filter gene eQTL file, {total_gene_file_count} gene files, good result {good_result_count}')
        return eqtl_filter_result

    def __prepare_eqtl_data(self, eqtl_trait_file_path, gene_id):
        try:
            eqtl_trait_df = pd.read_table(eqtl_trait_file_path, sep=const.column_spliter, header=0,
                                          usecols=self.qtl_col_dict.values(),
                                          dtype={self.qtl_col_dict['chrom']: 'category',
                                                 self.qtl_col_dict['position']: 'Int64'})
        except Exception as e:
            logging.error(f'error to prepare eqtl data: {eqtl_trait_file_path}')
            logging.error(f'exception: {e}')
        # eqtl_trait_df[Processor.VAR_ID_COL_NAME] = 'chr' + \
        #     eqtl_trait_df[self.qtl_col_dict['chrom']].astype(str) + \
        #     '_' + eqtl_trait_df[self.qtl_col_dict['position']].astype(str) + \
        #     '_' + eqtl_trait_df[self.qtl_col_dict['ref']].astype(str) + \
        #     '_' + eqtl_trait_df[self.qtl_col_dict['alt']].astype(str)

        if self.qtl_col_dict['variant_id'] not in set(eqtl_trait_df.columns):
            eqtl_trait_df[self.qtl_col_dict['variant_id']] = \
                'chr' + eqtl_trait_df[self.qtl_col_dict['chrom']].astype(str) + '_' + \
                eqtl_trait_df[self.qtl_col_dict['position']].astype(str)+ '_' + \
                eqtl_trait_df[self.qtl_col_dict['ref']].astype(str)+ '_' + \
                eqtl_trait_df[self.qtl_col_dict['alt']].astype(str)
                
        eqtl_trait_df[Processor.VAR_ID_COL_NAME] = eqtl_trait_df[self.qtl_col_dict['variant_id']].apply(lambda x: '_'.join(x.split('_')[:4]))
        utils.clean_data(eqtl_trait_df, dup_consider_subset=Processor.VAR_ID_COL_NAME, keep_dup='first')

        # utils.drop_indel_snp(eqtl_trait_df, self.qtl_col_dict['alt'], self.qtl_col_dict['ref'])

        if eqtl_trait_df.empty:
            del eqtl_trait_df
            utils.delete_file_if_exists(eqtl_trait_file_path)
            logging.info(f'No required data in eQTL file {eqtl_trait_file_path}')
            return None
        else:
            pval_filter_eqtl_df = \
                utils.filter_data_frame_by_p_value(eqtl_trait_df, self.config_holder.qtl_p_threshold,
                                                   self.qtl_col_dict['pvalue'], inplace=False)
            pval_filter_eqtl_df = pval_filter_eqtl_df.sort_values(self.qtl_col_dict['pvalue'])
            # pval_filter_eqtl_df = eqtl_trait_df.drop(
            #     eqtl_trait_df[(eqtl_trait_df[self.qtl_col_dict['pvalue']] == 0.0)].index,
            #     inplace=False)
            # logging.info(f"************qtl_p_threshold: {self.config_holder.qtl_p_threshold}")
            if pval_filter_eqtl_df.empty:
                del eqtl_trait_df
                del pval_filter_eqtl_df
                utils.delete_file_if_exists(eqtl_trait_file_path)
                return None
            else:
                chrom_list = []
                lead_SNP_list = []
                positions_list = []
                cluster_pos_dict = {}
                count = 0
                self.__qtl_ld_interval(pval_filter_eqtl_df, eqtl_trait_df, 
                                       self.qtl_window_size, 
                                       self.qtl_preprocesed_dir,
                                       self.qtl_col_dict, chrom_list, 
                                       lead_SNP_list, positions_list, 
                                       cluster_pos_dict, count, 
                                       self.config_holder.qtl_p_threshold,
                                       self.qtl_LD_window, gene_id,
                                       additional_expansion = self.qtl_LD_additional_expansion
                                       )
                eqtl_trait_df.sort_values([self.qtl_col_dict['chrom'], self.qtl_col_dict['position']],
                                          inplace=True)
                eqtl_trait_df.to_csv(eqtl_trait_file_path, sep=const.output_spliter, index=False)

                return eqtl_trait_df 



    def __qtl_ld_interval(self, pval_filter_df, total_df, window_size, vcf_output_path,
                          col_dict, chrom_list, lead_SNP_list, positions_list, 
                          cluster_pos_dict, count, p_threshold, outputpath, phenotype_id,
                          additional_expansion = 50000):

        population = self.global_config.get('population', 'EUR').upper()


        for _, row in pval_filter_df.iterrows(): 
            chrom = row.loc[col_dict['chrom']]
            pos = row.loc[col_dict['position']]
            chrom_cluster_pos_list = cluster_pos_dict.get(chrom, [])
            input_vcf = os.path.join(self.ref_vcf_dir, population, f'chr{chrom}.vcf.gz')
            # print(f"chrom {chrom}, pos {pos}, chrom_cluster_pos_list {chrom_cluster_pos_list}")
            if pos in set(chrom_cluster_pos_list):
                continue
            else:
                # chrom_cluster_pos_list.append(pos)
                # cluster_pos_dict[chrom] = chrom_cluster_pos_list
                # retrieve range_df from group by peak_positions.min <= group[position] <= peak_positions.max
                cluster_start = max(0, int(pos) - int(window_size))
                cluster_end = int(pos) + int(window_size)

                range_df = total_df[(total_df[col_dict['chrom']] == chrom) & (
                        cluster_start <= total_df[col_dict['position']]) & (
                                        total_df[col_dict['position']] <= cluster_end)]
                lead_SNP_id = \
                    range_df.loc[range_df[col_dict['position']] == pos][Processor.VAR_ID_COL_NAME].iloc[0]
                print(f"__qtl_ld_interval: {lead_SNP_id}") # chr1_13550_G_A

                vcf_output_dir = os.path.join(vcf_output_path, 'vcf')
                Path(vcf_output_dir).mkdir(parents=True, exist_ok=True)
                output_vcf_name = f"{lead_SNP_id}.vcf"

                utils.extract_vcf_data(chrom, range_df, input_vcf, vcf_output_dir,
                               output_vcf_name, col_dict['position'], 
                               target_snp_col_name=Processor.VAR_ID_COL_NAME, 
                               extract_step_size=window_size)
                output_ld_file = os.path.join(vcf_output_dir, f"{lead_SNP_id}")
                if lead_SNP_id == None:
                    continue
                else:
                    os.system(self.shell_command_plink_execute.format(f'{vcf_output_dir}/{output_vcf_name}', 
                                                                lead_SNP_id, window_size, self.qtl_LD_r2_filter,
                                                                output_ld_file))
                    
                if Path(f"{output_ld_file}.ld").exists():
                    ld_filter_df = pd.read_csv(f"{output_ld_file}.ld", sep='\s+')
                    if len(ld_filter_df) < 2: # No SNP with LD > threshold with the lead_SNP_id
                        logging.info(f"No SNP has a linkage disequilibrium (LD) with the target SNP {lead_SNP_id} that exceeds the {self.qtl_LD_r2_filter}.")
                        continue
                    else:
                        # ld_filter_df
                        range_start = int(list(ld_filter_df['BP_B'])[0])
                        range_end = int(list(ld_filter_df['BP_B'])[-1])

                        print(f"range: {range_start} - {range_end}")

                        range_df[col_dict['position']] = range_df[col_dict['position']].astype('int')
                        print(f"print(range_df[col_dict['position']].dtype): {range_df[col_dict['position']].dtype}")
                        range_df = range_df[
                            (range_start-additional_expansion <= range_df[col_dict['position']]) & 
                            (range_df[col_dict['position']] <= range_end+additional_expansion)]
                        

                        file_name = f'{lead_SNP_id}-chr{chrom}.tsv.gz'
                        if len(range_df[range_df[col_dict['pvalue']] < p_threshold]) < 2:
                            logging.info(f"Less than 2 significant variants in {lead_SNP_id} loci")
                            continue

                        chrom_list.append(chrom)
                        lead_SNP_list.append(lead_SNP_id)
                        print(f"chrom: {chrom}")
                        print(f"lead SNP id:{lead_SNP_id}")
                        count += 1
                        chrom_cluster_pos_list = chrom_cluster_pos_list + list(range_df[col_dict['position']])
                        cluster_pos_dict[chrom] = list(set(chrom_cluster_pos_list))
                        positions_list.append(range_df[col_dict['position']].tolist())
                        
                del range_df

        # print(f"length of chrom_list: {len(chrom_list)}, {chrom_list}")
        # print(f"length of lead_SNP_list: {len(lead_SNP_list)}, {lead_SNP_list}")
        # print(f"length of positions_list: {len(positions_list)}, {positions_list}")
        cluster_summary_df = pd.DataFrame(
            {'chrom': chrom_list, 'lead_SNP': lead_SNP_list, 'positions': positions_list})
        cluster_summary_df['phenotype_id'] = phenotype_id

        if os.path.exists(outputpath) and os.path.getsize(outputpath) > 0:
            mode = 'a'
            header = False
        else:
            mode = 'w'
            header = True

        # torus input file include columns variant_id、loc、zscore
        cluster_summary_df.to_csv(outputpath, sep=const.output_spliter, 
                                  mode=mode, header=header, index=False)




################################################################################
#                             1.3 Preprocess sQTL                              #
################################################################################
# phenotype_id                       pval_nominal  slope     slope_se     chr  pos       ref alt  maf       rs_id_dbSNP155_GRCh38p13  gene_id          var_id_
# chr1_52606945_52608262_clu_3971_+  0.28387       -0.0549   0.051209215  1    51603025  G   A    0.193749  rs4477330               ENSG00000116157.6  chr1_51603025_G_A
# chr1_52606945_52608262_clu_3971_+  0.97581       -0.0035   0.11825493   1    51603307  T   G    0.03125   rs17106697              ENSG00000116157.6  chr1_51603307_T_G

    def preprocess_sqtl(self):
        # group by trait
        logging.info('start to split sQTL file by phenotype')
        utils.delete_dir(self.qtl_grouped_dir)
        Path(self.qtl_grouped_dir).mkdir(exist_ok=True, parents=True)
        utils.split_file_by_col_name_sqtl(self.qtl_grouped_dir, 
                                     f'{self.global_config["input"]["qtl"]["file"]}',
                                     self.qtl_col_dict['chrom'], 
                                     self.qtl_col_dict['gene_id'], 
                                     self.qtl_col_dict['phenotype_id'], 
                                     readonly_cols=self.qtl_col_dict.values(), 
                                     sep=self.config_holder.qtl_sep)
        logging.info('finish split sQTL file')
        # drop na or no significant snp gene file, sort by position
        total_pheno_file_count = 0
        good_result_count = 0
        chrom_list = []
        gene_list = []
        pheno_file_list = []
        postions_list = []
        logging.info('start to filter gene sQTL file')
        logging.info(f'Filtering sQTL data by p-value threshold {self.config_holder.qtl_p_threshold}')
        for chrom_dir in os.listdir(f'{self.qtl_grouped_dir}'):
            for gene_dir in os.listdir(f'{os.path.join(self.qtl_grouped_dir, chrom_dir)}'):
                if not chrom_dir.startswith('.'):
                    for qtl_file in os.listdir(f'{self.qtl_grouped_dir}/{chrom_dir}/{gene_dir}'):
                        # if qtl_file.upper().startswith(const.gene_id_prefix):
                        total_pheno_file_count += 1
                        _, high_risk_snp_df = self.__prepare_sqtl_data(
                            f'{self.qtl_grouped_dir}/{chrom_dir}/{gene_dir}/{qtl_file}')
                        if high_risk_snp_df is not None and len(high_risk_snp_df) > 1:
                            # print(f"len(high_risk_snp_df): {len(high_risk_snp_df)}")
                            good_result_count += 1
                            chrom_list.append(chrom_dir)
                            gene_list.append(gene_dir)
                            pheno_file_list.append(qtl_file)
                            postions_list.append(high_risk_snp_df[self.qtl_col_dict['position']].tolist())
                            logging.debug(f'good gene {qtl_file}')
        sqtl_filter_result = pd.DataFrame(
            {'chrom': chrom_list, 'gene': gene_list, 'pheno_file': pheno_file_list, 'positions': postions_list})
        sqtl_filter_result.to_csv(self.qtl_output_report, sep=const.output_spliter, index=False)
        logging.info(
            f'finish filter gene sQTL file, {total_pheno_file_count} gene files, good result {good_result_count}')
        return sqtl_filter_result

    def __prepare_sqtl_data(self, sqtl_trait_file_path):
        try:
            sqtl_trait_df = pd.read_table(sqtl_trait_file_path, sep=const.column_spliter, header=0,
                                          usecols=self.qtl_col_dict.values(),
                                          dtype={self.qtl_col_dict['chrom']: 'category',
                                                 self.qtl_col_dict['position']: 'Int64'})
        except Exception as e:
            logging.error(f'error to prepare sqtl data: {sqtl_trait_file_path}')
            logging.error(f'exception: {e}')
        # sqtl_trait_df[Processor.VAR_ID_COL_NAME] = 'chr' + \
        #     sqtl_trait_df[self.qtl_col_dict['chrom']].astype(str) + \
        #     '_' + sqtl_trait_df[self.qtl_col_dict['position']].astype(str) + \
        #     '_' + sqtl_trait_df[self.qtl_col_dict['ref']].astype(str) + \
        #     '_' + sqtl_trait_df[self.qtl_col_dict['alt']].astype(str)


        
        sqtl_trait_df[Processor.VAR_ID_COL_NAME] = sqtl_trait_df[self.qtl_col_dict['variant_id']].apply(lambda x: '_'.join(x.split('_')[:4]))
        utils.clean_data(sqtl_trait_df, dup_consider_subset=Processor.VAR_ID_COL_NAME, keep_dup='first')
        # utils.drop_indel_snp(sqtl_trait_df, self.qtl_col_dict['alt'], self.qtl_col_dict['ref'])
        if sqtl_trait_df.empty:
            del sqtl_trait_df
            utils.delete_file_if_exists(sqtl_trait_file_path)
            logging.info(f'No required data in sQTL file {sqtl_trait_file_path}')
            return None, None
        else:
            pval_filter_sqtl_df = sqtl_trait_df.drop(
                sqtl_trait_df[(sqtl_trait_df[self.qtl_col_dict['pvalue']] == 0.0)].index,
                inplace=False)
            if pval_filter_sqtl_df.empty:
                del sqtl_trait_df
                del pval_filter_sqtl_df
                utils.delete_file_if_exists(sqtl_trait_file_path)
                return None, None
            else:
                sqtl_trait_df.sort_values([self.qtl_col_dict['chrom'], self.qtl_col_dict['position']],
                                          inplace=True)
                sqtl_trait_df.to_csv(sqtl_trait_file_path, sep=const.output_spliter, index=False)
                return sqtl_trait_df, pval_filter_sqtl_df






################################################################################
#                              1.4 Preprocess GWAS                             #
################################################################################

    def preprocess_gwas(self): # focus on target loci
        print(f"preprocess_gwas")
        utils.delete_dir(self.gwas_preprocessed_dir)
        Path(self.gwas_preprocessed_dir).mkdir(exist_ok=True, parents=True)
        gwas_file_path = self.global_config['input']['gwas']['file']
        population = self.global_config.get('population', 'EUR').upper()
        logging.info(f'Reading GWAS file {gwas_file_path}')

        # If the GWAS file have not been preprocessed
        if not Path(self.gwas_cluster_summary).exists(): 
            print(f"Starting preprocess GWAS ... ")
        ############################################################################
        #                              1.4.1 Clean GWAS                            #
        ############################################################################
            gwas_df = pd.read_csv(gwas_file_path, sep=self.config_holder.gwas_sep, header=0,
                                    # usecols=self.gwas_col_dict.values(),
                                    dtype={self.gwas_col_dict['position']: 'Int64', 
                                           self.gwas_col_dict['chrom']: 'category',
                                        # self.gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                        # self.gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)
                                        self.gwas_col_dict['effect_allele']: 'category',
                                        self.gwas_col_dict['other_allele']: 'category'
                                        })
            raw_gwas_snp_size = len(gwas_df)
            logging.info(f'GWAS data memory usage: {gwas_df.memory_usage(deep=True)}')
            gwas_df[self.gwas_col_dict['chrom']] = gwas_df[self.gwas_col_dict['chrom']].astype(str).str.lower().str.strip(
                'chr').astype('category')
            logging.info(f'GWAS data dropping non-autosome data, time: {datetime.datetime.now()}')
            gwas_df.drop(labels=gwas_df[~gwas_df[self.gwas_col_dict['chrom']].isin([str(i) for i in range(1, 23)])].index,
                        inplace=True)
            logging.info(f'GWAS data converting EA and NEA to upper case, time: {datetime.datetime.now()}')
            # gwas_df[self.gwas_col_dict['effect_allele']] = gwas_df[self.gwas_col_dict['effect_allele']].str.upper().astype(
            #     pd.CategoricalDtype(const.SNP_ALLELE)) # only ['A', 'T', 'C', 'G']
            # gwas_df[self.gwas_col_dict['other_allele']] = gwas_df[self.gwas_col_dict['other_allele']].str.upper().astype(
            #     pd.CategoricalDtype(const.SNP_ALLELE)) # only ['A', 'T', 'C', 'G']
            gwas_df[self.gwas_col_dict['effect_allele']] = gwas_df[self.gwas_col_dict['effect_allele']].str.upper()
            gwas_df[self.gwas_col_dict['other_allele']] = gwas_df[self.gwas_col_dict['other_allele']].str.upper()

            # logging.info(f'GWAS data dropping INDEL SNPs, time: {datetime.datetime.now()}')
            # utils.drop_indel_snp(gwas_df, self.gwas_col_dict['effect_allele'], self.gwas_col_dict['other_allele'])
            # logging.info(f'GWAS data cleaning, time: {datetime.datetime.now()}')
            # utils.clean_data(gwas_df, dup_consider_subset=[self.gwas_col_dict['chrom'], self.gwas_col_dict['position']])
            if self.gwas_col_dict['variant_id'] not in gwas_df.columns:
                # self.gwas_col_dict['variant_id'] = 'variant_id'
                gwas_df[self.gwas_col_dict['variant_id']] = \
                    'chr' + gwas_df[self.gwas_col_dict['chrom']].astype(str) + '_' + \
                    gwas_df[self.gwas_col_dict['position']].astype(str)+ '_' + \
                    gwas_df[self.gwas_col_dict['other_allele']].astype(str)+ '_' + \
                    gwas_df[self.gwas_col_dict['effect_allele']].astype(str) # chr1_13550 -> chr1_13550_G_A
            
            print(f"gwas variant_id {gwas_df[self.gwas_col_dict['variant_id']]}")
            utils.clean_data(gwas_df, dup_consider_subset=[self.gwas_col_dict['variant_id']], keep_dup='first')

            gwas_df.drop(index=gwas_df[gwas_df[self.gwas_col_dict['se']] == 0].index, inplace=True) # ？
            filtered_gwas_snp_size = len(gwas_df)
            logging.info(
                f'GWAS data total {filtered_gwas_snp_size} rows, sorting by chrom and pos, time: {datetime.datetime.now()}')
            gwas_df.sort_values([self.gwas_col_dict['chrom'], self.gwas_col_dict['position']], inplace=True)
            logging.info(f'Merging alt ref from vcf into GWAS data')
            # discontinuous index cost a lot more memory
            gwas_df.reset_index(drop=True, inplace=True)
            # Merge alt/ref from vcf into gwas file for later use
            gwas_df['alt_'] = gwas_df[self.gwas_col_dict['effect_allele']]
            gwas_df['ref_'] = gwas_df[self.gwas_col_dict['other_allele']]
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
            gwas_df.drop_duplicates(subset=[self.gwas_col_dict['variant_id']],
                                    keep=False, inplace=True)
            # Fill SNP ref/alt by other_allele/effect_allele if the SNP is not in vcf
            gwas_df['ref_'].mask(gwas_df['ref_'].isna(), gwas_df[self.gwas_col_dict['other_allele']], inplace=True)
            gwas_df['alt_'].mask(gwas_df['alt_'].isna(), gwas_df[self.gwas_col_dict['effect_allele']], inplace=True)

            gwas_df[Processor.VAR_ID_COL_NAME] = gwas_df[self.gwas_col_dict['variant_id']].apply(lambda x: '_'.join(x.split('_')[:4]))
            logging.info(
                f'Writing GWAS preprocessed data to {self.gwas_preprocessed_file}, time: {datetime.datetime.now()}')
            gwas_df.to_csv(self.gwas_preprocessed_file, sep=const.output_spliter, header=True, index=False)
            logging.info(f'Filtering GWAS data by p-value threshold {self.config_holder.gwas_p_threshold}')
        ############################################################################
        #                        1.4.2 Significant GWAS SNPs                       #
        ############################################################################
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
            lead_SNP_list = []
            positions_list = []
            cluster_pos_dict = {}
            count = 0
            gene_list = []

        ############################################################################
        #                    1.4.3 Three ways to preprocess GWAS                   #
        ############################################################################
            if self.gwas_genomic_window == 'LD_based_window':
                self.__preprocess_gwas_LD_based(pval_filter_gwas_df, gwas_df, 
                                                self.gwas_window_size, chrom_list, 
                                                lead_SNP_list, positions_list, 
                                                cluster_pos_dict, count, 
                                                raw_gwas_snp_size, self.target_loci, 
                                                filtered_gwas_snp_size, self.gwas_LD_additional_expansion)
            elif self.gwas_genomic_window == 'fixed_GWAS_Loci_window':
                self.__preprocess_gwas_loci_based(pval_filter_gwas_df, gwas_df, 
                                                  self.gwas_window_size, chrom_list, 
                                                  lead_SNP_list, positions_list, 
                                                  cluster_pos_dict, count, 
                                                  raw_gwas_snp_size, self.target_loci, 
                                                  filtered_gwas_snp_size)
            elif self.gwas_genomic_window == 'Gene_based_window':
                self.__preprocess_gwas_gene_based(pval_filter_gwas_df, gwas_df, 
                                                  self.gwas_window_size, chrom_list, 
                                                  lead_SNP_list, gene_list, positions_list, 
                                                  cluster_pos_dict, count, 
                                                  raw_gwas_snp_size, self.target_loci, 
                                                  filtered_gwas_snp_size)

            logging.info(f'GWAS data preprocess completed, time: {datetime.datetime.now()}')
            

        # If the GWAS file have  been preprocessed
        else:
            print(f"Preprocessed GWAS already exist ... ")
            if self.target_loci == 'ALL' or None:
                pass
            else:
                cluster_summary_df = pd.read_csv(self.gwas_cluster_summary, sep=const.output_spliter)
                target_loci_df = pd.read_csv(self.target_loci, sep=const.output_spliter)
                cluster_summary_df.index = cluster_summary_df['lead_SNP']
                cluster_summary_df = cluster_summary_df.loc[list(set(cluster_summary_df['lead_SNP'] & set(target_loci_df['lead_SNP'])))]
                cluster_summary_df = cluster_summary_df.sort_values('lead_SNP')
                cluster_summary_df.to_csv(self.gwas_cluster_summary, sep=const.output_spliter, header=True, index=False)

            logging.info(f'GWAS data preprocess completed, time: {datetime.datetime.now()}')

    ############################################################################
    #                          preprocess_gwas_LD_based                        #
    ############################################################################
    def __preprocess_gwas_LD_based(self, pval_filter_gwas_df, gwas_df, window_size, 
                                   chrom_list, lead_SNP_list, positions_list, 
                                   cluster_pos_dict, count, raw_gwas_snp_size, 
                                   target_loci, filtered_gwas_snp_size, additional_expansion = 50000):
        logging.info("Preprocess GWAS file based on LD ...")
        population = self.global_config.get('population', 'EUR').upper()

        # shell_command_plink_execute = 'plink --silent --vcf {} --r2 --matrix --mac 1 --write-snplist --out {}'

        shell_command_plink_execute = 'plink --silent --vcf {} --r2 --ld-snp {} --ld-window-kb {} --ld-window-r2 {} --out {}'
        for _, row in pval_filter_gwas_df.iterrows(): 
            chrom = row.loc[self.gwas_col_dict['chrom']]
            pos = row.loc[self.gwas_col_dict['position']]
            chrom_cluster_pos_list = cluster_pos_dict.get(chrom, [])
            input_vcf = os.path.join(self.ref_vcf_dir, population, f'chr{chrom}.vcf.gz')
            # for cp in chrom_cluster_pos_list:
            #     cp_start = max(0, cp - window_size)
            #     cp_end = cp + window_size
            #     if cp_start <= pos <= cp_end:
            #         break
            if pos in set(chrom_cluster_pos_list):
                continue
            else:
                # chrom_cluster_pos_list.append(pos)
                # cluster_pos_dict[chrom] = chrom_cluster_pos_list
                # retrieve range_df from group by peak_positions.min <= group[position] <= peak_positions.max
                cluster_start = max(0, pos - window_size)
                cluster_end = pos + window_size

                range_df = gwas_df[(gwas_df[self.gwas_col_dict['chrom']] == chrom) & (
                        cluster_start <= gwas_df[self.gwas_col_dict['position']]) & (
                                        gwas_df[self.gwas_col_dict['position']] <= cluster_end)]
                lead_SNP_id = \
                    range_df.loc[range_df[self.gwas_col_dict['position']] == pos][Processor.VAR_ID_COL_NAME].iloc[0]
                print(f"__preprocess_gwas_LD_based: {lead_SNP_id}")
                    # chr1_13550_G_A
                # target_snp = \
                #     range_df.loc[range_df[self.gwas_col_dict['position']] == pos][self.gwas_col_dict['snp']].iloc[0]
                vcf_output_dir = os.path.join(self.gwas_preprocessed_dir, 'vcf')
                Path(vcf_output_dir).mkdir(parents=True, exist_ok=True)
                output_vcf_name = f"{lead_SNP_id}.vcf"

                utils.extract_vcf_data(chrom, range_df, input_vcf, vcf_output_dir,
                               output_vcf_name, self.gwas_col_dict['position'], 
                               target_snp_col_name=Processor.VAR_ID_COL_NAME, 
                               extract_step_size=window_size)
                output_ld_file = os.path.join(vcf_output_dir, f"{lead_SNP_id}")
                if lead_SNP_id == None:
                    continue
                else:
                    os.system(self.shell_command_plink_execute.format(f'{vcf_output_dir}/{output_vcf_name}', 
                                                                lead_SNP_id, window_size, self.gwas_LD_r2_filter,
                                                                output_ld_file))
                    
                
                if not Path(f"{output_ld_file}.ld").exists():
                    continue

                ld_filter_df = pd.read_csv(f"{output_ld_file}.ld", sep='\s+')
                
                if len(ld_filter_df) == 0:
                    logging.info(f"Calculating LD failure {lead_SNP_id}.")
                    continue
                if len(ld_filter_df) == 1: # No SNP with LD > threshold with the lead_SNP_id
                    logging.info(f"No SNP has a linkage disequilibrium (LD) with the target SNP {lead_SNP_id} that exceeds the {self.gwas_LD_r2_filter}.")

                range_start = int(list(ld_filter_df['BP_B'])[0])
                range_end = int(list(ld_filter_df['BP_B'])[-1])

                logging.info(f"range: {range_start} - {range_end}")
                range_df[self.gwas_col_dict['position']] = range_df[self.gwas_col_dict['position']].astype('int')

                range_df = range_df[
                    (range_start-additional_expansion <= range_df[self.gwas_col_dict['position']]) & 
                    (range_df[self.gwas_col_dict['position']] <= range_end+additional_expansion)]

                if self.global_LD_based_window == True and Path(self.qtl_LD_window).exists():
                    # get the union of gwas and qtl LD based window
                    logging.info(f"Preprocessing global LD based window ... ")

                    ld_window = list(range_df[self.gwas_col_dict['position']])
                    
                    overlap = False

                    qtl_LD_window_summary = pd.read_csv(self.qtl_LD_window, sep=const.column_spliter)
                    for ix_qtl, row_qtl in qtl_LD_window_summary.iterrows():
                        chrom_qtl = str(row_qtl.loc['chrom'])
                        if str(chrom_qtl) != str(chrom):
                            continue

                        qtl_positions = ast.literal_eval(row_qtl.loc['positions'])
                        if len(set(range_df[self.gwas_col_dict['position']]) & set(qtl_positions)) < self.min_matching_number:
                            logging.info(f"Less than {self.min_matching_number} matched variants")
                            continue
                        else:
                            overlap = True
                            logging.info(f"More than {self.min_matching_number} matched variants")
                            ld_window = ld_window + list(set(qtl_positions))
                            logging.info(f"ld_window: {ld_window}")
                    
                    if overlap == False:
                        continue

                    range_df = gwas_df[(gwas_df[self.gwas_col_dict['chrom']] == chrom) & (
                            int(min(ld_window)) <= gwas_df[self.gwas_col_dict['position']]) & (
                                            gwas_df[self.gwas_col_dict['position']] <= int(max(ld_window)))]
                    logging.info(f"len of range_df: {range_df}")
                    file_name = f'{lead_SNP_id}-chr{chrom}.tsv.gz'
                    if len(range_df[range_df[self.gwas_col_dict['pvalue']] < self.config_holder.gwas_p_threshold]) < 2:
                        logging.info(f"Less than 2 significant variants in {lead_SNP_id} loci")
                        continue
                    logging.info(f"More than 2 significant variants in {lead_SNP_id} loci")
                    chrom_list.append(chrom)
                    lead_SNP_list.append(lead_SNP_id)
                    count += 1
                    range_df.to_csv(os.path.join(self.gwas_cluster_output_dir, file_name), 
                                    sep=const.output_spliter, header=True, index=False)
                    chrom_cluster_pos_list = chrom_cluster_pos_list + list(range_df[self.gwas_col_dict['position']])
                    cluster_pos_dict[chrom] = list(set(chrom_cluster_pos_list))
                    logging.info(f"cluster_pos_dict[chrom]: {cluster_pos_dict[chrom]}")
                    positions_list.append(range_df[self.gwas_col_dict['position']].tolist())

                else:
                    file_name = f'{lead_SNP_id}-chr{chrom}.tsv.gz'
                    if len(range_df[range_df[self.gwas_col_dict['pvalue']] < self.config_holder.gwas_p_threshold]) < 2:
                        logging.info(f"Less than 2 significant variants in {lead_SNP_id} loci")
                        continue
                    logging.info(f"More than 2 significant variants in {lead_SNP_id} loci")
                    chrom_list.append(chrom)
                    lead_SNP_list.append(lead_SNP_id)
                    count += 1
                    range_df.to_csv(os.path.join(self.gwas_cluster_output_dir, file_name), sep=const.output_spliter,
                                    header=True,
                                    index=False)
                    logging.info(f"range_dfrange_df: {range_df}")
                    chrom_cluster_pos_list = chrom_cluster_pos_list + list(range_df[self.gwas_col_dict['position']])
                    cluster_pos_dict[chrom] = list(set(chrom_cluster_pos_list))
                    logging.info(f"cluster_pos_dict[chrom]: {cluster_pos_dict[chrom]}")
                    positions_list.append(range_df[self.gwas_col_dict['position']].tolist())
                    del range_df
                #del pval_filter_range_df
        del pval_filter_gwas_df

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

        with open(f'{self.gwas_preprocessed_dir}/gwas_snp_summary.log', mode='w') as snp_file:
            snp_file.write(f'raw GWAS snp size: {raw_gwas_snp_size}\n')
            snp_file.write(f'filtered GWAS snp size: {filtered_gwas_snp_size}\n')
            snp_file.write(f'GWAS significant loci size: {count}\n')
        logging.info(
            f'Writing GWAS significant ranges summary to {self.gwas_cluster_summary}, time: {datetime.datetime.now()}')
        logging.info(f"chrom: {chrom_list}, lead_SNP: {lead_SNP_list}, positions: {positions_list}")
        cluster_summary_df = pd.DataFrame(
            {'chrom': chrom_list, 'lead_SNP': lead_SNP_list, 'positions': positions_list})
        

        ## User-defined target loci 
        if target_loci == 'ALL' or None:
            cluster_summary_df.to_csv(self.gwas_cluster_summary, sep=const.output_spliter, header=True, index=False)
        else:
            target_loci_df = pd.read_csv(target_loci, sep=const.output_spliter)
            cluster_summary_df.index = cluster_summary_df['lead_SNP']
            cluster_summary_df = cluster_summary_df.loc[list(set(cluster_summary_df['lead_SNP'] & set(target_loci_df['lead_SNP'])))]
            cluster_summary_df = cluster_summary_df.sort_values('lead_SNP')
            cluster_summary_df.to_csv(self.gwas_cluster_summary, sep=const.output_spliter, header=True, index=False)

        logging.info(f'GWAS data preprocess completed, time: {datetime.datetime.now()}')

        pass
    
    ############################################################################
    #                         preprocess_gwas_loci_based                       #
    ############################################################################
    def __preprocess_gwas_loci_based(self, pval_filter_gwas_df, gwas_df, window_size, 
                                   chrom_list, lead_SNP_list, positions_list, 
                                   cluster_pos_dict, count, raw_gwas_snp_size, 
                                   target_loci, filtered_gwas_snp_size):
        logging.info("Preprocess GWAS file based on GWAS loci ...")

        for _, row in pval_filter_gwas_df.iterrows(): 
            chrom = row.loc[self.gwas_col_dict['chrom']]
            pos = row.loc[self.gwas_col_dict['position']]
            chrom_cluster_pos_list = cluster_pos_dict.get(chrom, [])
            for cp in chrom_cluster_pos_list:
                cp_start = max(0, cp - window_size)
                cp_end = cp + window_size
                if cp_start <= pos <= cp_end:
                    break
            else:
                chrom_cluster_pos_list.append(pos)
                cluster_pos_dict[chrom] = chrom_cluster_pos_list
                # retrieve range_df from group by peak_positions.min <= group[position] <= peak_positions.max
                cluster_start = max(0, pos - window_size)
                cluster_end = pos + window_size
                range_df = gwas_df[(gwas_df[self.gwas_col_dict['chrom']] == chrom) & (
                        cluster_start <= gwas_df[self.gwas_col_dict['position']]) & (
                                        gwas_df[self.gwas_col_dict['position']] <= cluster_end)]
                lead_SNP_id = \
                    range_df.loc[range_df[self.gwas_col_dict['position']] == pos][Processor.VAR_ID_COL_NAME].iloc[0]
                    # chr1_13550_G_A
                print(f"__preprocess_gwas_loci_based: {lead_SNP_id}")
                chrom_list.append(chrom)
                lead_SNP_list.append(lead_SNP_id)
                count += 1
                file_name = f'{lead_SNP_id}-chr{chrom}.tsv.gz'
                range_df.to_csv(os.path.join(self.gwas_cluster_output_dir, file_name), sep=const.output_spliter,
                                header=True,
                                index=False)


                positions_list.append(range_df[self.gwas_col_dict['position']].tolist())
                del range_df
                #del pval_filter_range_df
        del pval_filter_gwas_df

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

        with open(f'{self.gwas_preprocessed_dir}/gwas_snp_summary.log', mode='w') as snp_file:
            snp_file.write(f'raw GWAS snp size: {raw_gwas_snp_size}\n')
            snp_file.write(f'filtered GWAS snp size: {filtered_gwas_snp_size}\n')
            snp_file.write(f'GWAS significant loci size: {count}\n')
        logging.info(
            f'Writing GWAS significant ranges summary to {self.gwas_cluster_summary}, time: {datetime.datetime.now()}')
        cluster_summary_df = pd.DataFrame(
            {'chrom': chrom_list, 'lead_SNP': lead_SNP_list, 'positions': positions_list})
        
        ## User-defined target loci 
        if target_loci == 'ALL' or None:
            cluster_summary_df.to_csv(self.gwas_cluster_summary, sep=const.output_spliter, header=True, index=False)
        else:
            target_loci_df = pd.read_csv(target_loci, sep=const.output_spliter)
            cluster_summary_df.index = cluster_summary_df['lead_SNP']
            cluster_summary_df = cluster_summary_df.loc[list(set(cluster_summary_df['lead_SNP'] & set(target_loci_df['lead_SNP'])))]
            cluster_summary_df = cluster_summary_df.sort_values('lead_SNP')
            cluster_summary_df.to_csv(self.gwas_cluster_summary, sep=const.output_spliter, header=True, index=False)

        logging.info(f'GWAS data preprocess completed, time: {datetime.datetime.now()}')
        pass
    
    def __preprocess_gwas_gene_based(self, pval_filter_gwas_df, gwas_df, window_size, 
                                   chrom_list, lead_SNP_list, gene_list, positions_list, 
                                   cluster_pos_dict, count, raw_gwas_snp_size, 
                                   target_loci, filtered_gwas_snp_size):
        qtl_process_path = self.qtl_preprocesed_dir
        qtl_process_path = self.qtl_grouped_dir
        qtl_chrom_ls = os.listdir(qtl_process_path)
        for chrom in qtl_chrom_ls:
            gwas_chrom = gwas_df[gwas_df[self.gwas_col_dict['chrom']] == str(chrom)]
            gwas_chrom.index = gwas_chrom[self.gwas_col_dict['snp']]
            if len(gwas_chrom) == 0:
                continue
            else:
                eqtl_chrom_gene_path = os.path.join(qtl_process_path, chrom)
                eqtl_gene_ls = os.listdir(eqtl_chrom_gene_path)
                for eqtl_gene_file in eqtl_gene_ls:
                    gene_name = eqtl_gene_file.rsplit('.', 2)[0]

                    eqtl_gene_df = pd.read_csv(os.path.join(eqtl_chrom_gene_path, eqtl_gene_file),sep='\t')
                    eqtl_gene_df.index = eqtl_gene_df['rsid']
                    if len(set(eqtl_gene_df.index) & set(gwas_chrom.index)) < 2:
                        continue
                    else:
                        range_df = gwas_chrom.loc[list(set(eqtl_gene_df.index) & set(gwas_chrom.index))]
                        range_df = range_df.sort_values(self.gwas_col_dict['pvalue'])
                        lead_SNP_id = str(range_df[Processor.VAR_ID_COL_NAME].iloc[0])
                        lead_SNP_list.append(lead_SNP_id)
                        # range_lead_var_id = range_df[Processor.VAR_ID_COL_NAME].iloc[0]
                        chrom_list.append(chrom)
                        gene_list.append(gene_name)
                        loci_count += 1
                        file_name = f'{gene_name}-chr{chrom}.tsv.gz'
                        range_df.to_csv(os.path.join(self.gwas_cluster_output_dir, file_name), sep=const.output_spliter,
                                        header=True,
                                        index=False)
                        positions_list.append(range_df[self.gwas_col_dict['position']].tolist())
                        del range_df

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

        with open(f'{self.gwas_preprocessed_dir}/gwas_snp_summary.log', mode='w') as snp_file:
            snp_file.write(f'raw GWAS snp size: {raw_gwas_snp_size}\n')
            snp_file.write(f'filtered GWAS snp size: {filtered_gwas_snp_size}\n')
            snp_file.write(f'GWAS significant loci size: {loci_count}\n')
        logging.info(
            f'Writing GWAS significant ranges summary to {self.gwas_cluster_summary}, time: {datetime.datetime.now()}')
        cluster_summary_df = pd.DataFrame(
            {'chrom': chrom_list, 'lead_variant': lead_SNP_list, 'gene_id': gene_list, 'positions': positions_list})
        cluster_summary_df.to_csv(self.gwas_cluster_summary, sep=const.output_spliter, header=True, index=False)
        logging.info(f'GWAS data preprocess completed, time: {datetime.datetime.now()}')


        pass

if __name__ == '__main__':
    processor = Processor()
    processor.preprocess_gwas()
    processor.preprocess_eqtl()
