import logging
import os
from datetime import datetime
from pathlib import Path
import sys
import pandas as pd
import yaml
from common import global_data_process as gdp, utils, constants as const
from fdr import prob_fdr

class Fastenloc:
    COLOC_TOOL_NAME = 'fastenloc'

    def __init__(self):
        logging.info('init Fastenloc')

    # def run(self, biological_context=None, working_dir=None,
    #         qtl_finemapping_file=None,
    #         qtl_output_report=None,
    #         output_torus_output_file=None,
    #         gwas_snp_count=None,
    #         tools_config_file=None):

    #     start_time = datetime.now()
    #     logging.info(f'fastenloc start at: {start_time}')

    #     shell_command_fastenloc_execute = 'fastenloc -e {} -g {} -t {} -tv {} -prefix {} {}'
    #     output_analyze_output_dir = f'{working_dir}/analyzed'
    #     fastenloc_params = utils.get_tools_params(self.COLOC_TOOL_NAME, tools_config_file=tools_config_file,
    #                                               param_prefix='-')

    #     Path(output_analyze_output_dir).mkdir(parents=True, exist_ok=True)
    #     com_str = shell_command_fastenloc_execute.format(qtl_finemapping_file,
    #                                                      f'{output_torus_output_file}.gz',
    #                                                      biological_context, gwas_snp_count,
    #                                                      f'{output_analyze_output_dir}/{biological_context}',
    #                                                      fastenloc_params)
    #     logging.info(f'fastenloc com_str: {com_str}')
    #     os.system(com_str)
    #     report_output_snp_tsv_file = self.__analyze_result(output_analyze_output_dir, biological_context, qtl_output_report, working_dir)

    #     logging.info(f'fastenloc complete at: {datetime.now()},duration: {datetime.now() - start_time}, with params {fastenloc_params}')
    #     return report_output_snp_tsv_file

    # def __analyze_result(self, output_dir, final_report_file, qtl_output_report, working_dir):
    #     # report_output_sig_file = f'{output_dir}/{final_report_file}.enloc.sig.out'
    #     report_output_sig_file = f'{output_dir}/{final_report_file}.enloc.gene.out'

    #     report_output_sig_tsv_file = f'{output_dir}/{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
    #     # sort sig file LCP(locus-level colocalization probability)
    #     # sort sig file GLCP(gene locus-level colocalization probability)
    #     # fdrthreshold_outfile = os.path.join(working_dir, 'analyzed', 'fdr_threshold.txt')
    #     if utils.file_exists(report_output_sig_file):
    #         df_output_sig = pd.read_csv(report_output_sig_file, sep='\s+')
    #         if len(df_output_sig) > 0:
    #             # df_output_sig.columns = ['Signal', 'Num_SNP', 'CPIP_qtl', 'CPIP_gwas_marginal',
    #             #                          'CPIP_gwas_qtl_prior',
    #             #                          'RCP', 'LCP']
    #             # df_output_sig.sort_values(by='LCP', ascending=False, inplace=True)
    #             # df_output_sig['gene_id'] = df_output_sig['Signal'].str.split(':').str[0]
    #             df_output_sig.columns = ['gene_id', 'GRCP', 'GLCP']
    #             df_output_sig.sort_values(by='GLCP', ascending=False, inplace=True)

    #             filtered_gene = pd.read_csv(qtl_output_report, sep=const.column_spliter)
    #             filtered_gene['gene_id'] = filtered_gene['pheno_file'].str.split('.').str[0]

    #             df_output_sig_mg = pd.merge(df_output_sig, filtered_gene[['gene_id','chrom']], 
    #                                         on=['gene_id'], how='left')
    #             df_output_sig_mg = df_output_sig_mg.round(4)
    #             df_output_sig_mg.to_csv(report_output_sig_tsv_file, sep=const.output_spliter, header=True, index=False)
                
                
    #             ## FDR threshold
    #             prob_thresh, notes = prob_fdr.calc_threshold_for_prob_rpt(report_output_sig_tsv_file, 'GLCP')
    #             # print(f"fastenlocthreshold: {prob_thresh}")
    #             # config = {
    #             #     'value': float(prob_thresh),
    #             #     'note': notes,
    #             # }
    #             # with open(fdrthreshold_outfile, 'w') as file:
    #             #     yaml.dump(config, file, default_flow_style=False, sort_keys=False)
    #         # else:
    #         #     ## FDR threshold
    #         #     config = {
    #         #         'value': 1,
    #         #         'note': "No result found",
    #         #     }
    #         #     with open(fdrthreshold_outfile, 'w') as file:
    #         #         yaml.dump(config, file, default_flow_style=False, sort_keys=False)

    #     # else:
    #     #     ## FDR threshold
    #     #     config = {
    #     #         'value': 1,
    #     #         'note': "No result found",
    #     #     }
    #     #     with open(fdrthreshold_outfile, 'w') as file:
    #     #         yaml.dump(config, file, default_flow_style=False, sort_keys=False)


    def run(self, biological_context=None, working_dir=None,

            qtl_output_report=None,
            summary_statistics_file = None,
            tools_config_file=None):
        start_time = datetime.now()
        logging.info(f'fastenloc start at: {start_time}')

        shell_command_fastenloc_execute = 'fastenloc -sum {} -prefix {}'
        output_analyze_output_dir = f'{working_dir}/analyzed'
        fastenloc_params = utils.get_tools_params(self.COLOC_TOOL_NAME, tools_config_file=tools_config_file,
                                                  param_prefix='-')

        Path(output_analyze_output_dir).mkdir(parents=True, exist_ok=True)
        com_str = shell_command_fastenloc_execute.format(f'{summary_statistics_file}',
                                                         f'{output_analyze_output_dir}/{biological_context}',
                                                         fastenloc_params)
        logging.info(f'fastenloc com_str: {com_str}')
        os.system(com_str)
        report_output_snp_tsv_file = self.__analyze_result(output_analyze_output_dir, biological_context, qtl_output_report, working_dir)

        logging.info(f'fastenloc complete at: {datetime.now()},duration: {datetime.now() - start_time}, with params {fastenloc_params}')
        return report_output_snp_tsv_file


    def __analyze_result(self, output_dir, biological_context, qtl_output_report, working_dir):
        # report_output_sig_file = f'{output_dir}/{final_report_file}.enloc.sig.out'
        report_output_sig_file = f'{output_dir}/{biological_context}.enloc.gene.out'
        report_output_sig_tsv_file = f'{output_dir}/{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        logging.info(f"fastenloc report: {report_output_sig_tsv_file}")
        
        if utils.file_exists(report_output_sig_file):
            df_output_sig = pd.read_csv(report_output_sig_file, sep='\s+')
            if len(df_output_sig) > 0:
                df_output_sig.columns = ['ix', 'GRCP', 'GLCP']
                df_output_sig['gene_id'] = df_output_sig['ix'].apply(lambda x: x.split('_')[0])
                df_output_sig['phenotype_id'] = df_output_sig['ix'].apply(lambda x: '_'.join(x.split('_')[1:7]))
                df_output_sig['lead_variant'] = df_output_sig['ix'].apply(lambda x: '_'.join(x.split('_')[7:11]))
                df_output_sig['chrom'] = df_output_sig['ix'].apply(lambda x: x.split('_')[1].strip('chr'))
                df_output_sig.sort_values(by='GRCP', ascending=False, inplace=True)
                # filtered_gene = pd.read_csv(qtl_output_report, sep=const.column_spliter)
                # filtered_gene['phenotype_id'] = filtered_gene['pheno_file'].str.split('.').str[0]
                # df_output_sig_mg = pd.merge(df_output_sig[['phenotype_id','lead_variant','GRCP','GLCP']], filtered_gene[['phenotype_id','chrom']], 
                #                             on=['phenotype_id'], how='left')
                # df_output_sig_mg = df_output_sig_mg.round(4)
                df_output_sig = df_output_sig[['chrom', 'gene_id', 'phenotype_id', 'lead_variant', 'GRCP', 'GLCP']]
                df_output_sig = df_output_sig.round(4)
                df_output_sig.to_csv(report_output_sig_tsv_file, sep=const.output_spliter, header=True, index=False)
                
        return report_output_sig_tsv_file


if __name__ == '__main__':
    fastenloc = Fastenloc()
    processor = gdp.Processor()
    _working_dir = os.path.join(processor.tool_parent_dir, fastenloc.COLOC_TOOL_NAME)
    fastenloc.run(biological_context=processor.biological_context,
                  working_dir=_working_dir,
                  qtl_finemapping_file=processor.global_config['input']['qtl_finemapping_file'],
                  qtl_output_report=processor.qtl_output_report)
