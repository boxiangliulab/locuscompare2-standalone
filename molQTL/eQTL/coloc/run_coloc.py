import ast
import concurrent
import json
import logging
import os
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from fdr import prob_fdr
import pandas as pd
import yaml

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

class Coloc:
    COLOC_TOOL_NAME = 'coloc'

    def __init__(self):
        logging.info('init eQTL Coloc')

    def run(self,
            working_dir=None,
            var_id_col_name=None,
            gwas_cluster_output_dir=None,
            gwas_cluster_summary=None,
            gwas_col_dict=None,
            gwas_sample_size=None,
            qtl_output_report=None,
            qtl_grouped_dir=None,
            qtl_col_dict=None,
            qtl_sample_size=None,
            gwas_type=None,
            qtl_type=None,
            parallel=False,
            tools_config=None,
            parallel_worker_num=10, 
            rank_dir = None,
            current_analysis_order = None, 
            total_numof_analyses = None, 
            whether_schedual = False,
            min_matching_number = 0,
            gwas_threshold = 5.0e-8,
            qtl_threshold = 5.0e-8):
        
        
        gwas_type_dict = {gwas_col_dict['chrom']: 'category',
                          gwas_col_dict['position']: 'Int64',
                        #   gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                        #   gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)
                          gwas_col_dict['effect_allele']: 'category',
                          gwas_col_dict['other_allele']: 'category'
                          }
        qtl_type_dict = {qtl_col_dict['chrom']: 'category',
                          qtl_col_dict['position']: 'Int64',
                        #   qtl_col_dict['alt']: pd.CategoricalDtype(const.SNP_ALLELE),
                        #   qtl_col_dict['ref']: pd.CategoricalDtype(const.SNP_ALLELE),
                          qtl_col_dict['alt']: 'category',
                          qtl_col_dict['ref']: 'category',
                          qtl_col_dict['phenotype_id']: 'category'
                          }
        start_time = datetime.now()
        Path(working_dir).mkdir(parents=True, exist_ok=True)
        logging.info(f'run_coloc start at: {start_time}')
        qtl_summary_df = pd.read_csv(qtl_output_report, sep=const.column_spliter,
                                      dtype={qtl_col_dict['chrom']: 'category'})
        gwas_summary_df = pd.read_csv(gwas_cluster_summary, sep=const.column_spliter,
                                      dtype={gwas_col_dict['chrom']: 'category'})
        gwas_cluster_snps_dict = self.__get_cluster_significant_snps_dict(gwas_summary_df)
        del gwas_summary_df
        # Put gwas range files in a list
        gwas_range_files = {}
        for gwas_range_file in os.listdir(gwas_cluster_output_dir):
            part_list = utils.get_file_name(gwas_range_file).split('-')
            if len(part_list) < 2 or 'chr' not in part_list[1]:
                continue
            chrom = part_list[1].strip('chr')
            range_files = gwas_range_files.get(chrom, [])
            range_files.append(os.path.join(gwas_cluster_output_dir, gwas_range_file))
            gwas_range_files[chrom] = range_files
        coloc_input_dir = os.path.join(working_dir, 'input')
        Path(coloc_input_dir).mkdir(parents=True, exist_ok=True)
        output_file = self.get_output_file(working_dir)
        Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
        Path(self.__get_output_dir(working_dir)).mkdir(parents=True, exist_ok=True)
        # Loop to process all QTL trait file
        gwas_chroms = gwas_range_files.keys()
        # logging.info(f"gwas_chroms: {gwas_chroms}")
        # logging.info(f"gwas_range_files: {gwas_range_files}")
        _p1, _p2, _p12 = self.__get_coloc_run_params(tools_config)
        # logging.info(f"_p1, _p2, _p12: {_p1}, {_p2}, {_p12}")

        total_len = len(qtl_summary_df)
        if parallel:
            with ThreadPoolExecutor(max_workers=parallel_worker_num) as executor:
                futures = []
                for ix, row in qtl_summary_df.iterrows():
                    if whether_schedual:
                        outputschedule(rownum=ix,
                            totalnum=total_len,
                            current_analysis_order = current_analysis_order,
                            total_numof_analyses=total_numof_analyses,
                            rank_dir=rank_dir)
                        
                    chrom = str(row.loc['chrom'])
                    print(f"chrom: {chrom}")
                    print(f"gwas_chroms: {gwas_chroms}")
                    if chrom not in gwas_chroms:
                        logging.info(f"{chrom} not in gwas_chroms")
                        print((f"{chrom} not in gwas_chroms"))
                        continue
                    qtl_gene_file = os.path.join(qtl_grouped_dir, chrom, row.loc['pheno_file'])
                    # logging.info(f"qtl_gene_file: {qtl_gene_file}")
                    phenotype_id = utils.get_eqtl_gene_name(qtl_gene_file)
                    # logging.info(f"phenotype_id: {phenotype_id}")
                    for gwas_range_file in gwas_range_files[chrom]:
                        range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]

                        pre_output_file = os.path.join(self.__get_output_dir(working_dir), f'{phenotype_id}_{range_lead_snp}.tsv')
                        if Path(pre_output_file).exists():
                            continue
                        logging.info(f"pass")
                        qtl_significant_positions = ast.literal_eval(row.loc['positions'])
                        
                        if range_lead_snp not in set(gwas_cluster_snps_dict.keys()):
                            logging.info(f"coloc range_lead_snp {range_lead_snp} not in set(gwas_cluster_snps_dict.keys())")
                            continue
                        if len(set(gwas_cluster_snps_dict[range_lead_snp]) & set(qtl_significant_positions)) < min_matching_number:
                            continue
                        logging.info(f"passa")

                        futures.append(executor.submit(self.process_gene, self.__get_output_dir(working_dir),
                                                    gwas_range_file, gwas_type_dict, gwas_col_dict, row,
                                                    qtl_gene_file, qtl_type_dict,
                                                    var_id_col_name, coloc_input_dir, phenotype_id, qtl_col_dict,
                                                    gwas_sample_size,
                                                    qtl_sample_size, gwas_type, qtl_type, _p1, _p2, _p12, min_matching_number, gwas_threshold, qtl_threshold))

                for future in concurrent.futures.as_completed(futures):
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error("".join(traceback.TracebackException.from_exception(exc).format()))

        else:
            for ix, row in qtl_summary_df.iterrows():
                if whether_schedual:
                    outputschedule(rownum=ix,
                        totalnum=total_len,
                        current_analysis_order = current_analysis_order,
                        total_numof_analyses=total_numof_analyses,
                        rank_dir=rank_dir)
                        
                chrom = str(row.loc['chrom'])
                if chrom not in gwas_chroms:
                    logging.info(f"{chrom} not in gwas_chroms")
                    continue
                qtl_gene_file = os.path.join(qtl_grouped_dir, chrom, row.loc['pheno_file'])
                # logging.info(f"qtl_gene_file: {qtl_gene_file}")
                phenotype_id = utils.get_eqtl_gene_name(qtl_gene_file)
                # logging.info(f"phenotype_id: {phenotype_id}")
                for gwas_range_file in gwas_range_files[chrom]:
                    range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
                    pre_output_file = os.path.join(self.__get_output_dir(working_dir), f'{phenotype_id}_{range_lead_snp}.tsv')
                    if Path(pre_output_file).exists():
                        continue
                    logging.info(f"pass")
                    qtl_significant_positions = ast.literal_eval(row.loc['positions'])
                    if range_lead_snp not in set(gwas_cluster_snps_dict.keys()):
                        logging.info(f"coloc range_lead_snp {range_lead_snp} not in set(gwas_cluster_snps_dict.keys())")
                        continue
                    if len(set(gwas_cluster_snps_dict[range_lead_snp]) & set(qtl_significant_positions)) < min_matching_number:
                        continue
                    logging.info(f"passa")
                    self.process_gene(self.__get_output_dir(working_dir), gwas_range_file, gwas_type_dict,
                                        gwas_col_dict, row, qtl_gene_file, qtl_type_dict,
                                        var_id_col_name, coloc_input_dir, phenotype_id, qtl_col_dict, gwas_sample_size,
                                        qtl_sample_size, gwas_type, qtl_type, _p1, _p2, _p12, min_matching_number, gwas_threshold, qtl_threshold)

        output_file = self.get_output_file(working_dir)
        self.__analyze_result(self.__get_output_dir(working_dir), output_file, working_dir)
        # fdrthreshold_outfile = os.path.join(working_dir, 'analyzed', 'fdr_threshold.txt')
        if not os.path.exists(output_file) or os.path.getsize(output_file) <= 0:
            # config = {
            #     'value': 1,
            #     'note': "No result found",
            # }
            # with open(fdrthreshold_outfile, 'w') as file:
            #     yaml.dump(config, file, default_flow_style=False, sort_keys=False)
            logging.warning(f'Process completed, duration {datetime.now() - start_time}, no result found')
        else:
            ## FDR threshold
            # prob_thresh, notes = prob_fdr.calc_threshold_for_prob_rpt(output_file, 'overall_H4')
            # config = {
            #     'value': float(prob_thresh),
            #     'note': notes,
            # }
            # with open(fdrthreshold_outfile, 'w') as file:
            #     yaml.dump(config, file, default_flow_style=False, sort_keys=False)
            logging.info(
                f'Process completed, duration {datetime.now() - start_time}, with params p1: {_p1} p2:{_p2} p12:{_p12}, check {output_file} for result!')
            
        return output_file

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

    def process_gene(self, output_dir, gwas_range_file, gwas_type_dict, gwas_col_dict, row, qtl_gene_file,
                     qtl_type_dict,
                     var_id_col_name, coloc_input_dir, phenotype_id, qtl_col_dict, gwas_sample_size,
                     qtl_sample_size, gwas_type, qtl_type, p1, p2, p12, min_matching_number, gwas_threshold, qtl_threshold):
        logging.info(f'COLOC process_gene')
        logging.info(f'gwas_threshold: {gwas_threshold}')
        logging.info(f'qtl_threshold: {qtl_threshold}')
        # print(f'COLOC process_gene')
        range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
        # print(f"eqtl coloc range_lead_snp: {range_lead_snp}")
        candidate_gwas_df = pd.read_table(gwas_range_file, sep=const.column_spliter,
                                          usecols=[
                                              var_id_col_name,
                                              gwas_col_dict['chrom'],
                                              gwas_col_dict['position'],
                                              gwas_col_dict['effect_allele'],
                                              gwas_col_dict['other_allele'],
                                              gwas_col_dict['beta'],
                                              gwas_col_dict['se'],
                                              gwas_col_dict['pvalue']],
                                          dtype=gwas_type_dict)
        if len(candidate_gwas_df) <= 1:
            logging.info(f'no candidate_gwas_df')
            return

        qtl_trait_df = pd.read_table(qtl_gene_file, sep=const.column_spliter,
                                      usecols=[
                                          var_id_col_name,
                                          qtl_col_dict['chrom'],
                                          qtl_col_dict['position'],
                                          qtl_col_dict['alt'],
                                          qtl_col_dict['ref'],
                                          qtl_col_dict['beta'],
                                          qtl_col_dict['se'],
                                          qtl_col_dict['pvalue'],
                                          qtl_col_dict['phenotype_id'],
                                          qtl_col_dict['maf']],
                                      dtype=qtl_type_dict)

        qtl_trait_df.drop(
            index=qtl_trait_df[~qtl_trait_df[var_id_col_name].isin(candidate_gwas_df[var_id_col_name])].index,
            inplace=True)
        # print(f"qtl_trait_df[var_id_col_name]: {qtl_trait_df[var_id_col_name]}")
        # print(f"candidate_gwas_df[var_id_col_name]: {candidate_gwas_df[var_id_col_name]}")
        
        if len(qtl_trait_df) < min_matching_number:
            logging.info(f'Skip, qtl no more than {min_matching_number}') # 有问题
            return
        # print(f"len(qtl_trait_df): {len(qtl_trait_df)}")
        utils.drop_non_intersect_rows(qtl_trait_df, var_id_col_name, candidate_gwas_df, var_id_col_name)
        if len(candidate_gwas_df) < min_matching_number:
            logging.info(f'Skip, gwas no more than {min_matching_number}')
            return
        coloc_gwas_input_path = os.path.join(coloc_input_dir, f'gwas_{phenotype_id}_{range_lead_snp}.tsv.gz')
        coloc_qtl_input_path = os.path.join(coloc_input_dir, f'qtl_{phenotype_id}_{range_lead_snp}.tsv.gz')
        output_file = os.path.join(output_dir, f'{phenotype_id}_{range_lead_snp}.tsv')
        print(f'coloc_gwas_input_path: {coloc_gwas_input_path}')
        print(f'coloc_qtl_input_path: {coloc_qtl_input_path}')
        print(f'output_file: {output_file}')

        # Adjust qtl beta sign according to gwas allele order
        candidate_gwas_df.reset_index(drop=True, inplace=True)
        qtl_trait_df.reset_index(drop=True, inplace=True)
        # utils.adjust_allele_order(candidate_gwas_df,
        #                           gwas_col_dict['effect_allele'],
        #                           gwas_col_dict['other_allele'],
        #                           gwas_col_dict['chrom'],
        #                           gwas_col_dict['position'],
        #                           qtl_trait_df,
        #                           ref_df_chrom_col_name=qtl_col_dict['chrom'],
        #                           ref_df_pos_col_name=qtl_col_dict['position'],
        #                           ref_df_alt_allele_col_name=qtl_col_dict['alt'],
        #                           ref_df_ref_allele_col_name=qtl_col_dict['ref'],
        #                           gbeta_col_name=gwas_col_dict['beta'])


        # Now gwas/qtl/vcf have the same num of rows, write candidate data to file
        # Reverse GWAS&qtl column mapping key-value and pass to R so that dataframe in R has fixed column names
        if ('varbeta' not in qtl_col_dict.keys() or qtl_col_dict.get('varbeta') is None) and (
                'se' in qtl_col_dict.keys() and qtl_col_dict.get('se') is not None):
            qtl_trait_df['varbeta'] = qtl_trait_df[qtl_col_dict['se']] ** 2
        qtl_trait_df.rename({v: k for k, v in qtl_col_dict.items()}, axis='columns', inplace=True)
        if ('varbeta' not in gwas_col_dict.keys() or gwas_col_dict.get('varbeta') is None) and (
                'se' in gwas_col_dict.keys() and gwas_col_dict.get('se') is not None):
            candidate_gwas_df['varbeta'] = candidate_gwas_df[gwas_col_dict['se']] ** 2
        candidate_gwas_df.rename({v: k for k, v in gwas_col_dict.items()}, axis='columns', inplace=True)



        if len(qtl_trait_df[qtl_trait_df['pvalue'] < qtl_threshold]) <= 0:
            logging.info(f'{phenotype_id}_{range_lead_snp} coloc no sig qtl_trait {qtl_threshold}')
            return
        if len(candidate_gwas_df[candidate_gwas_df['pvalue'] < gwas_threshold]) <= 0:
            logging.info(f'{phenotype_id}_{range_lead_snp} coloc no sig candidate_gwas {gwas_threshold}')
            return


        qtl_trait_df.to_csv(coloc_qtl_input_path, sep=const.output_spliter, header=True, index=False)
        del qtl_trait_df
        candidate_gwas_df.to_csv(coloc_gwas_input_path, sep=const.output_spliter, header=True, index=False)
        del candidate_gwas_df

        logging.info(f'Running coloc for significant SNP {range_lead_snp} on '
             f'GWAS input file {coloc_gwas_input_path} and qtl input file {coloc_qtl_input_path}')
        rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'coloc.R')
        # logging.info(f"{rscript_path}")
        logging.info(f'Rscript --no-save --no-restore {rscript_path} '
                     f'{output_file} {coloc_gwas_input_path} {coloc_qtl_input_path} {range_lead_snp} '
                     f'{gwas_sample_size} {qtl_sample_size} {gwas_type} {qtl_type} {p1} {p2} {p12}')
        os.system(f'Rscript --no-save --no-restore {rscript_path} '
                  f'{output_file} {coloc_gwas_input_path} {coloc_qtl_input_path} {range_lead_snp} '
                  f'{gwas_sample_size} {qtl_sample_size} {gwas_type} {qtl_type} {p1} {p2} {p12}')
        logging.info("{} completed".format(output_file))

        return output_file

    def __get_output_dir(self, working_dir):
        return os.path.join(working_dir, 'output')

    def __get_coloc_run_params(self, tools_config):
        params = utils.get_tools_params_dict(self.COLOC_TOOL_NAME, tools_config)
        _p1 = 1.0E-4 if params.get('p1') is None else params['p1']
        _p2 = 1.0E-4 if params.get('p2') is None else params['p2']
        _p12 = 1.0E-5 if params.get('p12') is None else params['p12']
        return _p1, _p2, _p12

    def __analyze_result(self, output_dir, final_result_file, working_dir):
        logging.info(f"output_dir: {output_dir}")
        logging.info(f"final_result_file: {final_result_file}")
        logging.info(f"working_dir: {working_dir}")
        single_result_list = []
        for single_result in os.listdir(output_dir):
            if not (single_result.endswith('.tsv') or single_result.endswith('.tsv.gz')):
                continue
            single_result_list.append(pd.read_table(os.path.join(output_dir, single_result)))
        if len(single_result_list) == 0:
            return
        report_df = pd.concat(single_result_list)
        report_df.sort_values(by=['overall_H4', 'SNP.PP.H4'], ascending=False, inplace=True)
        report_df.rename(columns={'phenotype_id': 'gene_id'}, inplace=True)
        report_df.drop_duplicates(subset=['snp', 'SNP.PP.H4', 'gene_id'], inplace=True)
        report_df = report_df.round(4)
        report_df.to_csv(final_result_file, sep=const.output_spliter, header=True, index=False)



    def get_output_file(self, working_dir):
        _output_file_name = f'{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        return os.path.join(working_dir, 'analyzed', _output_file_name)


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    if len(os.listdir(glob_processor.gwas_cluster_output_dir)) == 0 or not os.path.exists(
            glob_processor.qtl_output_report):
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, Coloc.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _qtl_sample_size = glob_processor.global_config['input']['qtl'].get('sample_size', 948)
    _gwas_type = glob_processor.global_config['input']['gwas'].get('type', 'cc')
    _qtl_type = glob_processor.global_config['input']['qtl'].get('type', 'quant')
    coloc = Coloc()
    coloc.run(_working_dir,
              gdp.Processor.VAR_ID_COL_NAME,
              glob_processor.gwas_cluster_output_dir,
              glob_processor.gwas_col_dict,
              _gwas_sample_size,
              glob_processor.qtl_output_report,
              glob_processor.qtl_grouped_dir,
              glob_processor.qtl_col_dict,
              _qtl_sample_size,
              _gwas_type,
              _qtl_type)
