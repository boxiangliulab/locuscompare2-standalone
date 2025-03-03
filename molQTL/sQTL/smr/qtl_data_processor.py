import ast
import concurrent
import logging
import os
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import utils, global_data_process as gdp, constants as const

def outputschedule(rownum, totalnum, current_analysis_order, total_numof_analyses, rank_dir):
    calculated_schedule = int(rownum/totalnum * 80/total_numof_analyses + 80/total_numof_analyses * (current_analysis_order - 1))
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()


class SmrQTLProcessor:
    # SMR_MIN_SNP = 2

    def __init__(self):
        logging.info('init SmrQTLProcessor')

    def prepare_ld_ref(self,
                       working_dir=None,
                       var_id_col_name=None,
                       gwas_filter_file=None,
                       gwas_col_dict=None,
                       qtl_output_report=None,
                       qtl_grouped_dir=None,
                       qtl_col_dict=None,
                       qtl_p_thresh=None,
                       ref_vcf_dir=None,
                       population=None,
                       rank_dir = None,
                       current_analysis_order = None,
                       total_numof_analyses = None,
                       whether_schedual = False,
                       parallel=False,
                       parallel_worker_num=2,
                       min_matching_number = 0):
        
        
        start_time = datetime.now()
        logging.info(f'Generating LD ref file at {start_time}')
        Path(self.__get_ldref_dir(working_dir)).mkdir(exist_ok=True, parents=True)
        Path(self.__get_output_vcf_dir(working_dir)).mkdir(exist_ok=True, parents=True)
        pval_filtered_gwas_df = pd.read_table(gwas_filter_file, sep=const.column_spliter,
                                              usecols=[gwas_col_dict['chrom'], gwas_col_dict['position']],
                                              dtype={gwas_col_dict['chrom']: 'category',
                                                     gwas_col_dict['position']: 'Int64'})
        qtl_summary_df = pd.read_csv(qtl_output_report, sep=const.column_spliter,
                                      dtype={qtl_col_dict['chrom']: 'category'})
        # Loop to process all sQTL trait file
        totalnumof_qtlloci = len(qtl_summary_df)
        gwas_chroms = pval_filtered_gwas_df[gwas_col_dict['chrom']].unique().tolist()
        if parallel:
            with ThreadPoolExecutor(max_workers=parallel_worker_num) as executor:
                futures = []
                for ix, row in qtl_summary_df.iterrows():
                    # print(f"row: {row}")
                    if whether_schedual == True:
                        outputschedule(rownum=ix,
                                    totalnum=totalnumof_qtlloci,
                                    current_analysis_order = current_analysis_order,
                                    total_numof_analyses=total_numof_analyses,
                                    rank_dir=rank_dir)
                    chrom = str(row.loc['chrom'])
                    # print(f"row: {row}")
                    if chrom not in gwas_chroms:
                        continue
                    qtl_gene_file = os.path.join(qtl_grouped_dir, chrom, row.loc['gene'], row.loc['pheno_file'])
                    # print(f"qtl_gene_file: {qtl_gene_file}")
                    pheno_id = utils.get_pheno_name(qtl_gene_file)
                    qtl_positions = ast.literal_eval(row.loc['positions'])
                    if len(qtl_positions) < min_matching_number:
                        continue
                    if not pval_filtered_gwas_df[pval_filtered_gwas_df[gwas_col_dict['chrom']] == chrom].loc[
                           :, gwas_col_dict['position']].isin(qtl_positions).any():
                        continue
                    input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chrom}.vcf.gz')
                    futures.append(executor.submit(self.gen_plink_binary_ld_ref, working_dir, qtl_gene_file,
                                                   var_id_col_name, qtl_col_dict, chrom, pheno_id, input_vcf,
                                                   qtl_p_thresh,
                                                   min_matching_number))
                for future in concurrent.futures.as_completed(futures):
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error("".join(traceback.TracebackException.from_exception(exc).format()))
        else:
            for ix, row in qtl_summary_df.iterrows():
                # print(f"row: {row}")
                if whether_schedual == True:
                    outputschedule(rownum=ix,
                                    totalnum=totalnumof_qtlloci,
                                    current_analysis_order = current_analysis_order,
                                    total_numof_analyses=total_numof_analyses,
                                    rank_dir=rank_dir)
                chrom = str(row.loc['chrom'])
                if chrom not in gwas_chroms:
                    continue
                qtl_gene_file = os.path.join(qtl_grouped_dir, chrom, row.loc['gene'], row.loc['pheno_file'])
                gene_id = row.loc['gene']
                pheno_id = utils.get_pheno_name(qtl_gene_file)
                qtl_positions = ast.literal_eval(row.loc['positions'])
                if len(qtl_positions) < min_matching_number:
                    continue
                if not pval_filtered_gwas_df[pval_filtered_gwas_df[gwas_col_dict['chrom']] == chrom].loc[
                       :, gwas_col_dict['position']].isin(qtl_positions).any():
                    continue
                input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chrom}.vcf.gz')
                self.gen_plink_binary_ld_ref(working_dir, qtl_gene_file, var_id_col_name, qtl_col_dict,
                                             chrom, pheno_id, gene_id, input_vcf, qtl_p_thresh,
                                             min_matching_number)
        logging.info(f'Generating LD ref file completed at {datetime.now()}, duration {datetime.now() - start_time}')
        return (
            self.__get_output_vcf_dir(working_dir), # processed/default/diag/Whole_Blood/EUR/sqtl/smr
            self.__get_output_vcf_file_pattern(), # chr{}_{}_{}.vcf
            self.__get_ldref_dir(working_dir), # processed/default/diag/Whole_Blood/EUR/sqtl/smr/ldref
            self.__get_ldref_file_pattern() # chr{}_{}_{}
        )

    def gen_plink_binary_ld_ref(self, working_dir, qtl_gene_file, var_id_col_name, qtl_col_dict,
                                chrom, pheno_id, gene_id, input_vcf, qtl_p_thresh,
                                min_matching_number):
        
        
        qtl_trait_df = pd.read_table(qtl_gene_file, sep=const.column_spliter,
                                      usecols=[var_id_col_name, qtl_col_dict['position'], qtl_col_dict['pvalue']],
                                      dtype={qtl_col_dict['position']: 'Int64'})
        # candidate_qtl_trait_df = utils.filter_data_frame_by_p_value(qtl_trait_df, qtl_p_thresh,
        #                                                              qtl_col_dict['pvalue'], inplace=False)
        # candidate_qtl_trait_df = qtl_trait_df.drop(
        #         qtl_trait_df[(qtl_trait_df[self.qtl_col_dict['pvalue']] == 0.0)].index,
        #         inplace=False)
        candidate_qtl_trait_df = qtl_trait_df
        del qtl_trait_df
        if candidate_qtl_trait_df.shape[0] == 0:
            logging.warning(f'{qtl_gene_file} does not contain significant records')
            return
        output_vcf_name = self.__get_output_vcf_file_pattern().format(chrom, gene_id, pheno_id)
        output_vcf_full_path = os.path.join(self.__get_output_vcf_dir(working_dir), output_vcf_name)

        if not os.path.exists(input_vcf):
            logging.warning(f'!ref vcf {input_vcf} does not exist')
            return
        logging.debug(f'Extracting vcf for gene {pheno_id}: {datetime.now()}')
        utils.extract_vcf_data(f'chr{chrom}', candidate_qtl_trait_df, input_vcf,
                               self.__get_output_vcf_dir(working_dir), output_vcf_name,
                               qtl_col_dict['position'], var_id_col_name)
        vcf_matching_file = os.path.join(self.__get_output_vcf_dir(working_dir), 'matching',
                                         f'chr{chrom}_{gene_id}_{pheno_id}.tsv')
        if not os.path.exists(vcf_matching_file) or os.path.getsize(vcf_matching_file) <= 0:
            logging.warning(f'No generated vcf file for gene {gene_id} {pheno_id}')
            return
        vcf_matching_df = pd.read_table(vcf_matching_file,
                                        sep=const.column_spliter, header=0,
                                        usecols=[qtl_col_dict['position'], var_id_col_name],
                                        dtype={qtl_col_dict['position']: 'Int64'})
        if len(vcf_matching_df) < min_matching_number:
            print(f"len(vcf_matching_df) < min_matching_number")
            return
        # Drop qtl rows that does not have vcf records
        utils.drop_non_intersect_rows(candidate_qtl_trait_df, var_id_col_name, vcf_matching_df,
                                      var_id_col_name)
        del vcf_matching_df
        if len(candidate_qtl_trait_df) < min_matching_number:
            print(f"len(candidate_qtl_trait_df) < min_matching_number")
            return
        del candidate_qtl_trait_df
        logging.debug(f'Generating plink binary LD ref for gene {gene_id} {pheno_id}: {datetime.now()}')
        # Generate ld ref file from subset vcf file
        output_ld_ref_path = os.path.join(self.__get_ldref_dir(working_dir),
                                          self.__get_ldref_file_pattern().format(chrom, gene_id, pheno_id))
        utils.delete_file_if_exists(f'{output_ld_ref_path}.bim')
        utils.delete_file_if_exists(f'{output_ld_ref_path}.fam')
        utils.delete_file_if_exists(f'{output_ld_ref_path}.bed')
        print(f'plink --silent --vcf {output_vcf_full_path} --make-bed --snps-only --out {output_ld_ref_path}')
        os.system(f'plink --silent --vcf {output_vcf_full_path} --make-bed --snps-only --out {output_ld_ref_path}')

    def __get_output_vcf_dir(self, working_dir):
        return os.path.join(working_dir, 'vcf')

    def __get_output_vcf_file_pattern(self):
        return 'chr{}_{}_{}.vcf'

    def __get_ldref_dir(self, working_dir):
        return os.path.join(working_dir, 'ldref')

    def __get_ldref_file_pattern(self):
        return 'chr{}_{}_{}'


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    if not os.path.exists(glob_processor.qtl_output_report) or os.path.getsize(
            glob_processor.qtl_output_report) <= 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process to preprocess qtl files?')
    if not os.path.exists(glob_processor.gwas_filter_file) or os.path.getsize(glob_processor.gwas_filter_file) <= 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process to preprocess gwas files?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, 'smr')
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    smr_qtl_processor = SmrQTLProcessor()
    smr_qtl_processor.prepare_ld_ref(_working_dir,
                                      gdp.Processor.VAR_ID_COL_NAME,
                                      glob_processor.gwas_filter_file,
                                      glob_processor.gwas_col_dict,
                                      glob_processor.qtl_output_report,
                                      glob_processor.qtl_grouped_dir,
                                      glob_processor.qtl_col_dict,
                                      glob_processor.config_holder.qtl_p_threshold,
                                      glob_processor.ref_vcf_dir,
                                      pop)
