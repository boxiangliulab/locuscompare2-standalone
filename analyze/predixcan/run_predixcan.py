import logging
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

import analyze.predixcan.gwas_data_processor as gpr

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import coloc_utils as utils, global_data_process as gdp, constants as const


class Predixcan:

    def __init__(self):
        logging.info('init Predixcan')

    def run(self,
            working_dir=None,
            model_db_path=None,
            prediction_snp_covariance_path=None,
            gwas_processed_file=None,
            gwas_col_dict=None,
            eqtl_output_report=None):
        start_time = datetime.now()
        output_file = self.__get_output_file(working_dir)
        Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
        predixan_path = shutil.which("SPrediXcan.py")
        if predixan_path is None:
            raise ValueError(f'SPrediXcan.py not found not found, '
                             f'did you config it in PATH env(or did you activate conda env?)')
        os.system(f'{sys.executable} {predixan_path} '
                  f'--model_db_path {model_db_path} '
                  f'--covariance {prediction_snp_covariance_path} '
                  f'--gwas_file {gwas_processed_file} '
                  f'--snp_column {gpr.PredixcanGwasProcessor.PREDIXCAN_VAR_ID_COL_NAME} '
                  f'--effect_allele_column {gwas_col_dict["effect_allele"]} '
                  f'--non_effect_allele_column {gwas_col_dict["other_allele"]} '
                  f'--beta_column {gwas_col_dict["beta"]} '
                  f'--se_column {gwas_col_dict["se"]} '
                  f'--overwrite '
                  f'--keep_non_rsid '
                  f'--model_db_snp_key varID '
                  f'--remove_ens_version '
                  f'--output_file {output_file}')
        if not os.path.exists(output_file) or os.path.getsize(output_file) <= 0:
            logging.warning(f'Process completed, duration {datetime.now() - start_time}, no result found')
            return output_file
        eqtl_summary_df = pd.read_table(eqtl_output_report, sep=const.column_spliter,
                                        usecols=['chrom', 'gene_file'],
                                        dtype={'chrom': 'string'})
        # gene column name should be the same as predixcan output gene column name
        eqtl_summary_df['gene'] = eqtl_summary_df['gene_file'].str.rstrip('.tsv.gz')
        eqtl_summary_df.drop(columns=[ 'gene_file'], inplace=True)
        # predixcan output column sep is comma
        result_df = pd.read_table(output_file, sep=',')
        result_df = pd.merge(left=result_df, right=eqtl_summary_df,
                             left_on='gene', right_on='gene',
                             how='left')
        result_df.rename(columns={'gene': 'gene_id'}, inplace=True)
        result_df.to_csv(output_file, sep=const.output_spliter, header=True, index=False)
        self.__analyze_result(output_file)
        if not os.path.exists(output_file) or os.path.getsize(output_file) <= 0:
            logging.warning(f'Process completed, duration {datetime.now() - start_time}, no result found')
        else:
            logging.info(
                f'Process completed, duration {datetime.now() - start_time}, check {output_file} for result!')
        return output_file

    def __get_predixcan_path(self, config_predixan_path):
        actual_predixcan_path = config_predixan_path
        if not config_predixan_path.endswith('SPrediXcan.py'):
            if config_predixan_path.endswith('software') or config_predixan_path.endswith('software/'):
                actual_predixcan_path = os.path.join(config_predixan_path, 'SPrediXcan.py')
            else:
                actual_predixcan_path = os.path.join(config_predixan_path, 'software', 'SPrediXcan.py')
        return actual_predixcan_path

    def __analyze_result(self, output_file_path):
        report_df = pd.read_table(output_file_path)
        report_df.dropna(subset='pvalue', inplace=True)
        report_df.sort_values(by='pvalue', ascending=True, inplace=True)
        if report_df.shape[0] == 0:
            utils.delete_file_if_exists(output_file_path)
        else:
            report_df.to_csv(output_file_path, sep=const.output_spliter, header=True, index=False)

    def __get_output_file(self, working_dir):
        _output_file_name = f'{gpr.PredixcanGwasProcessor.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        return os.path.join(working_dir, 'analyzed', _output_file_name)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise ValueError(f'The gwas processed file path argument is required!')
    _gwas_processed_file = sys.argv[1]
    if not os.path.exists(_gwas_processed_file) or os.path.getsize(_gwas_processed_file) <= 0:
        raise ValueError(f'Dependant files not found, did you run gwas_data_processor?')
    glob_processor = gdp.Processor()
    _working_dir = os.path.join(glob_processor.tool_parent_dir, gpr.PredixcanGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _model_db_path, _prediction_snp_covariance_path = utils.get_predixcan_ref_files(
        glob_processor.config_holder.global_config)
    predixcan = Predixcan()
    predixcan.run(_working_dir,
                  _model_db_path,
                  _prediction_snp_covariance_path,
                  _gwas_processed_file,
                  glob_processor.gwas_col_dict,
                  glob_processor.eqtl_output_report)
