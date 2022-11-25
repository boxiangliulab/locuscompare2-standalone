import ast
import os
import sys
import analyze.smr.gwas_data_processor as gpr
from pathlib import Path
from datetime import datetime
import pandas as pd
import logging
sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import coloc_utils as utils, global_data_process as gdp, constants as const


class SmrEqtlProcessor:
    SMR_MIN_SNP = 2

    def __init__(self):
        logging.info('init SmrEqtlProcessor')

    def prepare_ld_ref(self,
                       working_dir=None,
                       var_id_col_name=None,
                       gwas_filter_file=None,
                       gwas_col_dict=None,
                       eqtl_output_report=None,
                       eqtl_output_dir=None,
                       eqtl_col_dict=None,
                       eqtl_p_thresh=None,
                       ref_vcf_dir=None,
                       population=None):
        start_time = datetime.now()
        logging.info(f'Generating LD ref file at {start_time}')
        Path(self.__get_ldref_dir(working_dir)).mkdir(exist_ok=True, parents=True)
        Path(self.__get_output_vcf_dir(working_dir)).mkdir(exist_ok=True, parents=True)
        pval_filtered_gwas_df = pd.read_table(gwas_filter_file, sep=const.column_spliter,
                                              dtype={gwas_col_dict['chrom']: 'category',
                                                     gwas_col_dict['position']: 'Int64'})
        eqtl_summary_df = pd.read_csv(eqtl_output_report, sep=const.column_spliter,
                                      dtype={eqtl_col_dict['chrom']: 'category'})
        # Loop to process all eQTL trait file
        gwas_chroms = pval_filtered_gwas_df[gwas_col_dict['chrom']].unique().tolist()
        for _, row in eqtl_summary_df.iterrows():
            chrom = str(row.loc['chrom'])
            if chrom not in gwas_chroms:
                continue
            loop_start_time = datetime.now()
            eqtl_gene_file = os.path.join(eqtl_output_dir, chrom, row.loc['gene_file'])
            gene_id = utils.get_file_name(eqtl_gene_file)
            eqtl_positions = ast.literal_eval(row.loc['positions'])
            if len(eqtl_positions) < SmrEqtlProcessor.SMR_MIN_SNP:
                continue
            if not pval_filtered_gwas_df[pval_filtered_gwas_df[gwas_col_dict['chrom']] == chrom].loc[
                   :, gwas_col_dict['position']].isin(eqtl_positions).any():
                continue
            eqtl_trait_df = pd.read_table(eqtl_gene_file, sep=const.column_spliter,
                                          dtype={eqtl_col_dict['chrom']: 'category',
                                                 eqtl_col_dict['position']: 'Int64'})
            candidate_eqtl_trait_df = utils.filter_data_frame_by_p_value(eqtl_trait_df, eqtl_p_thresh,
                                                                         eqtl_col_dict['pvalue'], inplace=False)
            del eqtl_trait_df
            output_vcf_name = self.__get_output_vcf_file_pattern().format(chrom, gene_id)
            output_vcf_full_path = os.path.join(self.__get_output_vcf_dir(working_dir), output_vcf_name)
            input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chrom}.vcf.gz')
            if not os.path.exists(input_vcf):
                logging.warning(f'!ref vcf {input_vcf} does not exist')
                continue
            logging.debug(f'Extracting vcf for gene {gene_id}: {datetime.now()}')
            utils.extract_vcf_data(f'chr{chrom}', candidate_eqtl_trait_df, input_vcf,
                                   self.__get_output_vcf_dir(working_dir), output_vcf_name,
                                   eqtl_col_dict['position'], var_id_col_name)
            vcf_matching_file = os.path.join(self.__get_output_vcf_dir(working_dir), 'matching',
                                             f'chr{chrom}_{gene_id}.tsv')
            if not os.path.exists(vcf_matching_file) or os.path.getsize(vcf_matching_file) <= 0:
                logging.warning(f'No generated vcf file for gene {gene_id}')
                continue
            vcf_matching_df = pd.read_table(vcf_matching_file,
                                            sep=const.column_spliter, header=0,
                                            usecols=[eqtl_col_dict['position'], var_id_col_name],
                                            dtype={eqtl_col_dict['position']: 'Int64'})
            if len(vcf_matching_df) < SmrEqtlProcessor.SMR_MIN_SNP:
                continue
            # Drop eQTL rows that does not have vcf records
            utils.drop_non_intersect_rows(candidate_eqtl_trait_df, var_id_col_name, vcf_matching_df,
                                          var_id_col_name)
            del vcf_matching_df
            if len(candidate_eqtl_trait_df) < SmrEqtlProcessor.SMR_MIN_SNP:
                continue
            logging.debug(f'Generating plink binary LD ref for gene {gene_id}: {datetime.now()}')
            # Generate ld ref file from subset vcf file
            output_ld_ref_path = os.path.join(self.__get_ldref_dir(working_dir),
                                              self.__get_ldref_file_pattern().format(chrom, gene_id))
            utils.delete_file_if_exists(f'{output_ld_ref_path}.bim')
            utils.delete_file_if_exists(f'{output_ld_ref_path}.fam')
            utils.delete_file_if_exists(f'{output_ld_ref_path}.bed')
            os.system(f'plink -vcf {output_vcf_full_path} --make-bed --snps-only --out {output_ld_ref_path}')
            logging.debug(f'Generating LD ref file for gene {gene_id} done, duration {datetime.now() - loop_start_time}')
        logging.info(f'Generating LD ref file completed at {datetime.now()}, duration {datetime.now() - start_time}')
        return (
            self.__get_output_vcf_dir(working_dir),
            self.__get_output_vcf_file_pattern(),
            self.__get_ldref_dir(working_dir),
            self.__get_ldref_file_pattern()
        )

    def __get_output_vcf_dir(self, working_dir):
        return os.path.join(working_dir, 'vcf')

    def __get_output_vcf_file_pattern(self):
        return 'chr{}_{}.vcf'

    def __get_ldref_dir(self, working_dir):
        return os.path.join(working_dir, 'ldref')

    def __get_ldref_file_pattern(self):
        return 'chr{}_{}'


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    g_processor = gpr.SmrGwasProcessor()
    if not os.path.exists(glob_processor.eqtl_output_report) or os.path.getsize(
            glob_processor.eqtl_output_report) <= 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process to preprocess eqtl files?')
    if not os.path.exists(glob_processor.gwas_filter_file) or os.path.getsize(glob_processor.gwas_filter_file) <= 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process to preprocess gwas files?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, gpr.SmrGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _eqtl_p_thresh = glob_processor.global_config['p-value_threshold']['eqtl']
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    smr_eqtl_processor = SmrEqtlProcessor()
    smr_eqtl_processor.prepare_ld_ref(_working_dir,
                                      gdp.Processor.VAR_ID_COL_NAME,
                                      glob_processor.gwas_filter_file,
                                      glob_processor.gwas_col_dict,
                                      glob_processor.eqtl_output_report,
                                      glob_processor.eqtl_output_dir,
                                      glob_processor.eqtl_col_dict,
                                      _eqtl_p_thresh,
                                      glob_processor.ref_vcf_dir,
                                      pop)
