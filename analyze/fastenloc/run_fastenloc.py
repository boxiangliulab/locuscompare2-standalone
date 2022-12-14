import pandas as pd
from common import global_data_process as gdp, coloc_utils as utils, constants as const
import os
from datetime import datetime
from pathlib import Path
import analyze.fastenloc.gwas_data_processor as gdr
import logging


class Fastenloc:
    COLOC_TOOL_NAME = 'fastenloc'

    def __init__(self):
        logging.info('init Fastenloc')

    def run(self, eqtl_tissue=None, working_dir=None,
            eqtl_finemapping_file=None, gwas_preprocessed_file=None,
            eqtl_output_dir=None, var_id_col_name=None,
            gwas_col_dict=None, eqtl_output_report=None,
            eqtl_col_dict=None,output_torus_output_file=None):
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
        report_output_snp_tsv_file = self.__analyze_result(output_analyze_output_dir, eqtl_tissue, gwas_preprocessed_file, eqtl_output_dir,
                              var_id_col_name, gwas_col_dict, eqtl_output_report, eqtl_col_dict)

        logging.info(f'fastenloc complete at: {datetime.now()},duration: {datetime.now() - start_time}')
        return report_output_snp_tsv_file

    def __analyze_result(self, output_dir, final_report_file, gwas_preprocessed_file, eqtl_output_dir, var_id_col_name,
                         gwas_col_dict, eqtl_output_report, eqtl_col_dict):
        report_output_sig_file = f'{output_dir}/{final_report_file}.enloc.sig.out'
        report_output_snp_file = f'{output_dir}/{final_report_file}.enloc.snp.out'

        report_output_sig_tsv_file = f'{output_dir}/{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        report_output_snp_tsv_file = f'{output_dir}/output_{datetime.now().strftime("%Y%m%d%H%M%S")}.sig.tsv.gz'

        # sort sig file rcp(colocalization probability)
        if utils.file_exists(report_output_sig_file):
            df_output_sig = pd.read_csv(report_output_sig_file, sep='\s+')
            df_output_sig.columns = ['Signal', 'Num_SNP', 'CPIP_qtl', 'CPIP_gwas_marginal',
                                     'CPIP_gwas_qtl_prior',
                                     'RCP','LCP']
            df_output_sig.sort_values(by='RCP', ascending=False, inplace=True)
            df_output_sig['gene_id'] = df_output_sig['Signal'].str.split(':').str[0]

            filtered_gene = pd.read_csv(eqtl_output_report, sep=const.column_spliter)
            filtered_gene['gene_id'] = filtered_gene['gene_file'].str.split('.').str[0]

            df_output_sig_mg = pd.merge(df_output_sig, filtered_gene[['gene_id', 'chrom']], on=['gene_id'], how='left')

            df_output_sig_mg['eqtl_path'] = eqtl_output_dir + '/' + df_output_sig_mg['chrom'].astype(str) + '/' + \
                                            df_output_sig_mg['Signal'].str.split(':').str[0] + '.tsv.gz'
            df_output_sig_mg['gwas_path'] = gwas_preprocessed_file
            df_output_sig_mg['rsid'] = ''
            df_output_sig_mg.to_csv(report_output_sig_tsv_file, sep=const.output_spliter, header=True, index=False)

            # sort snp file rcp(colocalization probability)
            df_output_snp = pd.read_csv(report_output_snp_file, sep='\s+')
            df_output_snp.columns = ['Signal', 'SNP', 'PIP_qtl', 'PIP_gwas_marginal',
                                     'PIP_gwas_qtl_prior',
                                     'SCP']
            df_output_snp.sort_values(by='SCP', ascending=False, inplace=True)

            df_output_snp['chrom'] = df_output_snp['SNP'].str.split('_').str[0].str.replace('chr', '')
            df_output_snp[var_id_col_name] = df_output_snp['SNP'].str.split('_').str[0] + '_' + \
                                             df_output_snp['SNP'].str.split('_').str[1]
            df_output_snp['eqtl_path'] = eqtl_output_dir + '/' + df_output_snp['chrom'] + '/' + \
                                         df_output_snp['Signal'].str.split(':').str[0] + '.tsv.gz'
            df_output_snp['gwas_path'] = gwas_preprocessed_file
            df_output_snp['gene_id'] = df_output_snp['Signal'].str.split(':').str[0]
            # preprocessed_pd = pd.read_csv(gwas_preprocessed_file, sep='\s+')

            merge_pd = utils.mapping_var_id_to_rsid(df_output_snp, var_id_col_name, 'gene_id',
                                                    gwas_preprocessed_file, var_id_col_name,
                                                    gwas_col_dict, eqtl_col_dict)

            # merge_pd = pd.merge(df_output_snp, preprocessed_pd[[var_id_col_name, gwas_col_dict['snp']]],
            #                     on=[var_id_col_name],
            #                     how='left')

            merge_pd.to_csv(report_output_snp_tsv_file, sep=const.output_spliter, header=True, index=False)
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
