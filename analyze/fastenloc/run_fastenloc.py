import pandas as pd
from common import global_data_process as gdp, coloc_utils as utils, constants as const
import os
from datetime import datetime
from pathlib import Path
import logging


class Fastenloc:
    COLOC_TOOL_NAME = 'fastenloc'

    def __init__(self):
        logging.info('init Fastenloc')

    def run(self, eqtl_tissue=None, working_dir=None,
            eqtl_finemapping_file=None, gwas_preprocessed_file=None,
            eqtl_output_dir=None, var_id_col_name=None,
            gwas_col_dict=None, eqtl_output_report=None,
            eqtl_col_dict=None, output_torus_output_file=None):
        start_time = datetime.now()
        logging.info(f'fastenloc start at: {start_time}')

        shell_command_fastenloc_execute = 'fastenloc -e {} -g {} -t {} -prefix {} {}'
        output_analyze_output_dir = f'{working_dir}/analyzed'
        fastenloc_params = utils.get_tools_params(self.COLOC_TOOL_NAME, param_prefix='-')

        Path(output_analyze_output_dir).mkdir(parents=True, exist_ok=True)
        com_str = shell_command_fastenloc_execute.format(eqtl_finemapping_file,
                                                         f'{output_torus_output_file}.gz',
                                                         eqtl_tissue,
                                                         f'{output_analyze_output_dir}/{eqtl_tissue}',
                                                         fastenloc_params)
        os.system(com_str)
        report_output_snp_tsv_file = self.__analyze_result(output_analyze_output_dir, eqtl_tissue,
                                                           gwas_preprocessed_file, eqtl_output_dir,
                                                           var_id_col_name, gwas_col_dict, eqtl_output_report,
                                                           eqtl_col_dict)

        logging.info(f'fastenloc complete at: {datetime.now()},duration: {datetime.now() - start_time}')
        return report_output_snp_tsv_file

    def __analyze_result(self, output_dir, final_report_file, gwas_preprocessed_file, eqtl_output_dir, var_id_col_name,
                         gwas_col_dict, eqtl_output_report, eqtl_col_dict):
        report_output_sig_file = f'{output_dir}/{final_report_file}.enloc.sig.out'

        report_output_sig_tsv_file = f'{output_dir}/{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        # sort sig file rcp(colocalization probability)
        if utils.file_exists(report_output_sig_file):
            df_output_sig = pd.read_csv(report_output_sig_file, sep='\s+')
            if len(df_output_sig) > 0:
                df_output_sig.columns = ['Signal', 'Num_SNP', 'CPIP_qtl', 'CPIP_gwas_marginal',
                                         'CPIP_gwas_qtl_prior',
                                         'RCP', 'LCP']
                df_output_sig.sort_values(by='RCP', ascending=False, inplace=True)
                df_output_sig['gene_id'] = df_output_sig['Signal'].str.split(':').str[0]

                filtered_gene = pd.read_csv(eqtl_output_report, sep=const.column_spliter)
                filtered_gene['gene_id'] = filtered_gene['gene_file'].str.split('.').str[0]

                df_output_sig_mg = pd.merge(df_output_sig, filtered_gene[['gene_id', 'chrom']], on=['gene_id'],
                                            how='left')
                df_output_sig_mg.to_csv(report_output_sig_tsv_file, sep=const.output_spliter, header=True, index=False)
        return report_output_sig_tsv_file


if __name__ == '__main__':
    fastenloc = Fastenloc()
    processor = gdp.Processor()
    _working_dir = os.path.join(processor.tool_parent_dir, fastenloc.COLOC_TOOL_NAME)
    fastenloc.run(eqtl_tissue=processor.eqtl_tissue,
                  working_dir=_working_dir,
                  eqtl_finemapping_file=processor.global_config['input']['eqtl_finemapping_file'],
                  gwas_preprocessed_file=processor.gwas_preprocessed_file,
                  eqtl_output_dir=processor.eqtl_output_dir,
                  var_id_col_name=processor.VAR_ID_COL_NAME,
                  gwas_col_dict=processor.gwas_col_dict,
                  eqtl_output_report=processor.eqtl_output_report,
                  eqtl_col_dict=processor.eqtl_col_dict)