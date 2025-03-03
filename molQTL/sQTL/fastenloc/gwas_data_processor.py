import logging
import os
from datetime import datetime
from pathlib import Path
import sys
import pandas as pd
from common import global_data_process as gdp, constants as const, utils
import ast
import concurrent
import json
import logging
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import yaml


def outputschedule(rownum, totalnum, current_analysis_order, total_numof_analyses, rank_dir):
    calculated_schedule = int(rownum/totalnum * 80/total_numof_analyses + 80/total_numof_analyses * (current_analysis_order - 1))
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()


# # new version fastenloc & full gene
# class FastenlocGwasProcessor:

#     def __init__(self):
#         logging.info('init FastenlocGwasProcessor')

#     def __loc_array_init(self, chr_str, chr_len):
#         arr = [f'loc_{chr_str}_{x}' for x in range(chr_len)]
#         return arr

#     def get_output_gwas_dir(self, output_base_dir):
#         return f'{output_base_dir}/gwas'

#     def get_output_torus_output(self, output_gwas_dir):
#         return f'{output_gwas_dir}/pip'

#     def get_output_torus_output_file(self, output_torus_output_dir):
#         return f'{output_torus_output_dir}/torus_output.pip'

#     def prepare_fastenlocinput_data(self, working_dir=None, 
#                                     gwas_preprocessed_file=None,
#                                     gwas_col_dict=None, eqtl_file = None, 
#                                     qtl_col_dict = None, rank_dir=None,
#                                     current_analysis_order=None, total_numof_analyses=None, 
#                                     whether_schedual=False):
#         start_time = datetime.now()
#         logging.info(f'Grouping fastenloc preprocessed gwas file at {start_time}')

#         output_base_dir = working_dir
#         output_gwas_dir = self.get_output_gwas_dir(output_base_dir)
#         output_prefastenloc_dir = f'{output_gwas_dir}/fastenloc_input'
#         output_prefastenloc_file = f'{output_prefastenloc_dir}/pre_fastenloc_summarystat.tsv.gz'
#         utils.delete_file_if_exists(output_prefastenloc_file)


#         gwas_use_col = [gwas_col_dict['snp'], gwas_col_dict['beta'],
#                        gwas_col_dict['se']]
        

#         eqtl_use_col = [qtl_col_dict['chrom'], qtl_col_dict['phenotype_id'], qtl_col_dict['snp'], qtl_col_dict['beta'], qtl_col_dict['se']]
        
#         eqtl_df = pd.read_csv(eqtl_file, sep=const.column_spliter, usecols=eqtl_use_col)

#         Path(output_prefastenloc_dir).mkdir(parents=True, exist_ok=True)


#         gwas_preprocessed_df = pd.read_table(gwas_preprocessed_file, sep=const.column_spliter, usecols=gwas_use_col)


#         eqtl_df.index = eqtl_df[qtl_col_dict['se']]
#         gwas_preprocessed_df.index = gwas_preprocessed_df[gwas_col_dict['snp']]



#         # # 获取共同的 SNP
#         common_snp = list(set(gwas_preprocessed_df.index) & set(eqtl_df.index))
#         # # 筛选共同的 SNP
#         eqtl_df = eqtl_df.loc[common_snp]
#         gwas_preprocessed_df = gwas_preprocessed_df.loc[common_snp]
#         # # 合并 eqtl 和 gwas 数据

#         merged_df = pd.merge(eqtl_df[[qtl_col_dict['phenotype_id'], qtl_col_dict['snp'], qtl_col_dict['beta'], qtl_col_dict['se']]], 
#                              gwas_preprocessed_df[[gwas_col_dict['beta'], gwas_col_dict['se']]], 
#                              left_index=True, 
#                              right_index=True, 
#                              how='left')

#         # 保存合并的数据到指定路径
#         merged_df.to_csv(output_prefastenloc_file, sep=const.output_spliter, header=False, index=False, compression='gzip')

        
#         if whether_schedual == True:
#             outputschedule(current_analysis_order = current_analysis_order,
#                         total_numof_analyses=total_numof_analyses,
#                         rank_dir=rank_dir)

#         # os.system(f"cp {output_torus_output_file}.gz ~/scratch/torus.output")
#         logging.info(f'prepare fastenloc gwas file at: {datetime.now()}, duration: {datetime.now() - start_time}')
#         return output_prefastenloc_file



# new version fastenloc & gwas loci
class FastenlocGwasProcessor:
    def __init__(self):
        logging.info('init FastenlocGwasProcessor')

    def __loc_array_init(self, chr_str, chr_len):
        arr = [f'loc_{chr_str}_{x}' for x in range(chr_len)]
        return arr

    def get_output_preprocess_dir(self, output_base_dir):
        return f'{output_base_dir}/fastenloc_input'


    def __convert_positions_str_to_list(self, positions_str):

        if isinstance(positions_str, str):
            return json.loads(positions_str)
        elif isinstance(positions_str, list):
            return positions_str
        else:
            return []

    def __get_cluster_significant_snps_dict(self, cluster_df):

        cluster_snps_dict = {}
        for _, row in cluster_df.iterrows():
            cluster_snps_dict[row.loc['lead_SNP']] = self.__convert_positions_str_to_list(row.loc['positions'])
        return cluster_snps_dict
    
    def prepare_fastenlocinput_data(self, 
                                    working_dir=None,
                                    var_id_col_name=None,
                                    gwas_cluster_output_dir=None,
                                    gwas_cluster_summary=None,
                                    gwas_col_dict=None,
                                    qtl_output_report=None,
                                    qtl_grouped_dir=None,
                                    qtl_col_dict=None,
                                    parallel=False,
                                    parallel_worker_num=10, 
                                    rank_dir = None,
                                    current_analysis_order = None, 
                                    total_numof_analyses = None, 
                                    whether_schedual = False,
                                    min_matching_number = 0):
        


        start_time = datetime.now()
        Path(working_dir).mkdir(parents=True, exist_ok=True)
        logging.info(f'Grouping fastenloc preprocessed gwas file at {start_time}')

        qtl_summary_df = pd.read_csv(qtl_output_report, sep=const.column_spliter,
                                      dtype={qtl_col_dict['chrom']: 'category'})
        gwas_summary_df = pd.read_csv(gwas_cluster_summary, sep=const.column_spliter,
                                      dtype={gwas_col_dict['chrom']: 'category'})
        gwas_cluster_snps_dict = self.__get_cluster_significant_snps_dict(gwas_summary_df)
        del gwas_summary_df
        # Put gwas range files in a list
        gwas_range_files = {}
        for gwas_range_file in os.listdir(gwas_cluster_output_dir):
            # print(f"gwas_range_file: {gwas_range_file}")
            part_list = utils.get_file_name(gwas_range_file).split('-')
            if len(part_list) < 2 or 'chr' not in part_list[1]:
                continue
            chrom = part_list[1].strip('chr')
            range_files = gwas_range_files.get(chrom, [])
            range_files.append(os.path.join(gwas_cluster_output_dir, gwas_range_file))
            gwas_range_files[chrom] = range_files

        # Loop to process all eQTL trait file
        gwas_chroms = gwas_range_files.keys()

        total_len = len(qtl_summary_df)

        output_base_dir = working_dir
        output_preprocess_dir = self.get_output_preprocess_dir(output_base_dir)
        output_prefastenloc_file = f'{output_preprocess_dir}/pre_fastenloc_summarystat.tsv.gz'
        output_prefastenloc_mapping_file = f'{output_preprocess_dir}/pre_fastenloc_mapping.tsv.gz'
        Path(output_preprocess_dir).mkdir(parents=True, exist_ok=True)
        utils.delete_file_if_exists(output_prefastenloc_file)
        utils.delete_file_if_exists(output_prefastenloc_mapping_file)

        for ix, row in qtl_summary_df.iterrows():
            # print(f"row: {row}")
            if whether_schedual:
                outputschedule(rownum=ix,
                    totalnum=total_len,
                    current_analysis_order = current_analysis_order,
                    total_numof_analyses=total_numof_analyses,
                    rank_dir=rank_dir)
                    
            chrom = str(row.loc['chrom'])
            if chrom not in gwas_chroms:
                continue
            sqtl_pheno_file = os.path.join(qtl_grouped_dir, chrom, row.loc['gene'], row.loc['pheno_file'])
            pheno_id = utils.get_pheno_name(sqtl_pheno_file)
            gene_id = row.loc['gene']
            for gwas_range_file in gwas_range_files[chrom]:
                eqtl_significant_positions = ast.literal_eval(row.loc['positions'])
                range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
                if len(set(gwas_cluster_snps_dict[range_lead_snp]) & set(eqtl_significant_positions)) == 0:
                    continue
                print(f"fffastenlocpro")
                self.process_pheno(working_dir, gwas_range_file,
                                    gwas_col_dict, row, sqtl_pheno_file,
                                    var_id_col_name, qtl_col_dict, 
                                    output_prefastenloc_file, 
                                    output_prefastenloc_mapping_file, 
                                    pheno_id, gene_id, min_matching_number)



        logging.info(f'prepare fastenloc gwas file at: {datetime.now()}, duration: {datetime.now() - start_time}')
        return output_prefastenloc_file



    def process_pheno(self, working_dir, gwas_range_file, gwas_col_dict, row, 
                      sqtl_pheno_file,var_id_col_name, qtl_col_dict, 
                      output_prefastenloc_file, output_prefastenloc_mapping_file, 
                      pheno_id, gene_id, min_matching_number):
        print(f"fastenloc sqtl process pheno")
        range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
        # print(f"range_lead_snp: {range_lead_snp}")
        candidate_gwas_df = pd.read_table(gwas_range_file, sep=const.column_spliter)
        if len(candidate_gwas_df) <= min_matching_number:
            print('fastenloc gwas fail')
            return
        qtl_trait_df = pd.read_table(sqtl_pheno_file, sep=const.column_spliter)
        # print(f"1: {qtl_trait_df[var_id_col_name][:5]}")
        qtl_trait_df.drop(
            index=qtl_trait_df[~qtl_trait_df[var_id_col_name].isin(candidate_gwas_df[var_id_col_name])].index,
            inplace=True)
        
        # print(f"2: {candidate_gwas_df[var_id_col_name][:5]}")
        if len(qtl_trait_df) <= min_matching_number:
            print('fastenloc no overlap')
            return
        utils.drop_non_intersect_rows(qtl_trait_df, var_id_col_name, candidate_gwas_df, var_id_col_name)
        if len(candidate_gwas_df) <= 1:
            return


        # candidate_gwas_df.reset_index(drop=True, inplace=True)
        # qtl_trait_df.reset_index(drop=True, inplace=True)

        candidate_gwas_df.index = candidate_gwas_df[var_id_col_name]
        # print(f"candidate_gwas_df index: {candidate_gwas_df.index}")
        qtl_trait_df.index = qtl_trait_df[var_id_col_name]
        qtl_trait_df['ix'] = f'{gene_id}_{pheno_id}_{range_lead_snp}'

        # merged_df = pd.merge(qtl_trait_df[['ix', qtl_col_dict['snp'], qtl_col_dict['beta'], qtl_col_dict['se']]], 
        #                      candidate_gwas_df[[gwas_col_dict['beta'], gwas_col_dict['se']]], 
        #                      left_index=True, 
        #                      right_index=True, 
        #                      how='left')
        merged_df = pd.merge(qtl_trait_df[['ix', qtl_col_dict['variant_id'], qtl_col_dict['beta'], qtl_col_dict['se'], qtl_col_dict['snp']]], 
                             candidate_gwas_df[[gwas_col_dict['beta'], gwas_col_dict['se']]], 
                             left_index=True, 
                             right_index=True, 
                             how='left')
        print(f"fastenloc merged_df: {merged_df.head()}")
        merged_df.columns = ['ix','variant_id','qtlbeta','qtlse','rsid','gwasbeta','gwasse']

        mapping_df = merged_df[['ix', 'variant_id','rsid']]
        merged_df = merged_df[['ix','variant_id','qtlbeta','qtlse','gwasbeta','gwasse']]

        if len(merged_df) < min_matching_number:
            logging.info(f"Warning: Skip Gene {gene_id} {pheno_id} {range_lead_snp} loci, no more than {min_matching_number} SNPs")
            return
        
        if os.path.exists(output_prefastenloc_file) and os.path.getsize(output_prefastenloc_file) > 0:
            mode = 'a'
            header = False
        else:
            mode = 'w'
            header = False

        # torus input file include columns variant_id、loc、zscore
        merged_df.to_csv(output_prefastenloc_file, sep=const.output_spliter, 
                         mode=mode, header=header, index=False)
        mapping_df.to_csv(output_prefastenloc_mapping_file, sep=const.output_spliter, 
                          mode=mode, header=header, index=False)
        print(f"output_prefastenloc_mapping_file: {output_prefastenloc_mapping_file}")


if __name__ == '__main__':
    fastenloc_gwas_pro = FastenlocGwasProcessor()
    processor = gdp.Processor()
    _working_dir = os.path.join(processor.tool_parent_dir, 'fastenloc')
    fastenloc_gwas_pro.prepare_gwas_data(working_dir=_working_dir,
                                         gwas_preprocessed_file=processor.gwas_preprocessed_file,
                                         gwas_col_dict=processor.gwas_col_dict,
                                         ld_block_loci_file=processor.global_config['input']['ld_block_loci_file'])
