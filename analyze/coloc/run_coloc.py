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

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import coloc_utils as utils, global_data_process as gdp, constants as const

def outputschedule(rownum, totalnum, currenttissuenum, numoftissues, rank_dir):
    calculated_schedule = int(rownum/totalnum * 80/numoftissues + 80/numoftissues * (currenttissuenum - 1))
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
        logging.info('init Coloc')

    def run(self,
            working_dir=None,
            var_id_col_name=None,
            gwas_cluster_output_dir=None,
            gwas_cluster_summary=None,
            gwas_col_dict=None,
            gwas_sample_size=None,
            eqtl_output_report=None,
            eqtl_output_dir=None,
            eqtl_col_dict=None,
            eqtl_sample_size=None,
            gwas_type=None,
            eqtl_type=None,
            parallel=False,
            tools_config=None,
            parallel_worker_num=10, 
            rank_dir = None,
            currenttissuenum = None, 
            numoftissues = None, 
            whether_schedual = False):
        
        
        gwas_type_dict = {gwas_col_dict['chrom']: 'category',
                          gwas_col_dict['position']: 'Int64',
                          gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                          gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)
                          }
        eqtl_type_dict = {eqtl_col_dict['chrom']: 'category',
                          eqtl_col_dict['position']: 'Int64',
                          eqtl_col_dict['alt']: pd.CategoricalDtype(const.SNP_ALLELE),
                          eqtl_col_dict['ref']: pd.CategoricalDtype(const.SNP_ALLELE),
                          eqtl_col_dict['gene_id']: 'category'
                          }
        start_time = datetime.now()
        Path(working_dir).mkdir(parents=True, exist_ok=True)
        logging.info(f'run_coloc start at: {start_time}')
        eqtl_summary_df = pd.read_csv(eqtl_output_report, sep=const.column_spliter,
                                      dtype={eqtl_col_dict['chrom']: 'category'})
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
        # Loop to process all eQTL trait file
        gwas_chroms = gwas_range_files.keys()
        _p1, _p2, _p12 = self.__get_coloc_run_params(tools_config)

        total_len = len(eqtl_summary_df)
        if parallel:
            with ThreadPoolExecutor(max_workers=parallel_worker_num) as executor:
                futures = []
                for ix, row in eqtl_summary_df.iterrows():
                    if whether_schedual:
                        outputschedule(rownum=ix,
                            totalnum=total_len,
                            currenttissuenum = currenttissuenum,
                            numoftissues=numoftissues,
                            rank_dir=rank_dir)
                        
                    chrom = str(row.loc['chrom'])
                    if chrom not in gwas_chroms:
                        continue
                    eqtl_gene_file = os.path.join(eqtl_output_dir, chrom, row.loc['gene_file'])
                    gene_id = utils.get_file_name(eqtl_gene_file)
                    for gwas_range_file in gwas_range_files[chrom]:
                        eqtl_significant_positions = ast.literal_eval(row.loc['positions'])
                        range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
                        if len(set(gwas_cluster_snps_dict[range_lead_snp]) & set(eqtl_significant_positions)) == 0:
                            continue
                        futures.append(executor.submit(self.process_gene, self.__get_output_dir(working_dir),
                                                       gwas_range_file, gwas_type_dict, gwas_col_dict, row,
                                                       eqtl_gene_file, eqtl_type_dict,
                                                       var_id_col_name, coloc_input_dir, gene_id, eqtl_col_dict,
                                                       gwas_sample_size,
                                                       eqtl_sample_size, gwas_type, eqtl_type, _p1, _p2, _p12))

                for future in concurrent.futures.as_completed(futures):
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error("".join(traceback.TracebackException.from_exception(exc).format()))

        else:
            for ix, row in eqtl_summary_df.iterrows():
                if whether_schedual:
                    outputschedule(rownum=ix,
                        totalnum=total_len,
                        currenttissuenum = currenttissuenum,
                        numoftissues=numoftissues,
                        rank_dir=rank_dir)
                        
                chrom = str(row.loc['chrom'])
                if chrom not in gwas_chroms:
                    continue
                eqtl_gene_file = os.path.join(eqtl_output_dir, chrom, row.loc['gene_file'])
                gene_id = utils.get_file_name(eqtl_gene_file)
                for gwas_range_file in gwas_range_files[chrom]:
                    eqtl_significant_positions = ast.literal_eval(row.loc['positions'])
                    range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
                    if len(set(gwas_cluster_snps_dict[range_lead_snp]) & set(eqtl_significant_positions)) == 0:
                        continue
                    self.process_gene(self.__get_output_dir(working_dir), gwas_range_file, gwas_type_dict,
                                      gwas_col_dict, row, eqtl_gene_file, eqtl_type_dict,
                                      var_id_col_name, coloc_input_dir, gene_id, eqtl_col_dict, gwas_sample_size,
                                      eqtl_sample_size, gwas_type, eqtl_type, _p1, _p2, _p12)

        self.__analyze_result(self.__get_output_dir(working_dir), output_file, working_dir)
        if not os.path.exists(output_file) or os.path.getsize(output_file) <= 0:
            logging.warning(f'Process completed, duration {datetime.now() - start_time}, no result found')
        else:
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
            cluster_snps_dict[row.loc['range_lead']] = self.__convert_positions_str_to_list(row.loc['positions'])
        return cluster_snps_dict

    def process_gene(self, output_dir, gwas_range_file, gwas_type_dict, gwas_col_dict, row, eqtl_gene_file,
                     eqtl_type_dict,
                     var_id_col_name, coloc_input_dir, gene_id, eqtl_col_dict, gwas_sample_size,
                     eqtl_sample_size, gwas_type, eqtl_type, p1, p2, p12):

        range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
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
            return
        eqtl_trait_df = pd.read_table(eqtl_gene_file, sep=const.column_spliter,
                                      usecols=[
                                          var_id_col_name,
                                          eqtl_col_dict['chrom'],
                                          eqtl_col_dict['position'],
                                          eqtl_col_dict['alt'],
                                          eqtl_col_dict['ref'],
                                          eqtl_col_dict['beta'],
                                          eqtl_col_dict['se'],
                                          eqtl_col_dict['pvalue'],
                                          eqtl_col_dict['gene_id'],
                                          eqtl_col_dict['maf']],
                                      dtype=eqtl_type_dict)
        eqtl_trait_df.drop(
            index=eqtl_trait_df[~eqtl_trait_df[var_id_col_name].isin(candidate_gwas_df[var_id_col_name])].index,
            inplace=True)
        if len(eqtl_trait_df) <= 1:
            return
        utils.drop_non_intersect_rows(eqtl_trait_df, var_id_col_name, candidate_gwas_df, var_id_col_name)
        if len(candidate_gwas_df) <= 1:
            return
        coloc_gwas_input_path = os.path.join(coloc_input_dir, f'gwas_{gene_id}_{range_lead_snp}.tsv.gz')
        coloc_eqtl_input_path = os.path.join(coloc_input_dir, f'eqtl_{gene_id}_{range_lead_snp}.tsv.gz')
        output_file = os.path.join(output_dir, f'{gene_id}_{range_lead_snp}.tsv')
        # output_vcf_path = os.path.join(vcf_output_dir, f'{range_lead_snp}-chr{chrom}.vcf')
        # vcf_matching_file = os.path.join(vcf_output_dir, 'matching',
        #                                  f'{range_lead_snp}-chr{chrom}.tsv')
        # if not os.path.exists(vcf_matching_file) or os.path.getsize(vcf_matching_file) <= 0:
        #     logging.info(f'No generated vcf file for significant SNP {range_lead_snp}')
        #     continue
        # vcf_matching_df = pd.read_table(vcf_matching_file,
        #                                 sep=const.column_spliter, header=0,
        #                                 usecols=[gwas_col_dict['position'], var_id_col_name],
        #                                 dtype={gwas_col_dict['position']: 'Int64'})
        # if len(vcf_matching_df) <= 1:
        #     continue
        # # Drop GWAS rows that does not have vcf records
        # # SNP in vcf_matching_df is subset of SNP in candidate_gwas_df, so it's fine to drop intersect rows here
        # utils.drop_non_intersect_rows(candidate_gwas_df, var_id_col_name, vcf_matching_df,
        #                               var_id_col_name)
        # # Drop eQTL rows that does not have vcf records
        # utils.drop_non_intersect_rows(eqtl_trait_df, var_id_col_name, vcf_matching_df,
        #                               var_id_col_name)
        # del vcf_matching_df
        # if len(candidate_gwas_df) <= 1:
        #     continue
        # Adjust eqtl beta sign according to gwas allele order
        candidate_gwas_df.reset_index(drop=True, inplace=True)
        eqtl_trait_df.reset_index(drop=True, inplace=True)
        utils.adjust_allele_order(candidate_gwas_df,
                                  gwas_col_dict['effect_allele'],
                                  gwas_col_dict['other_allele'],
                                  gwas_col_dict['chrom'],
                                  gwas_col_dict['position'],
                                  eqtl_trait_df,
                                  ref_df_chrom_col_name=eqtl_col_dict['chrom'],
                                  ref_df_pos_col_name=eqtl_col_dict['position'],
                                  ref_df_alt_allele_col_name=eqtl_col_dict['alt'],
                                  ref_df_ref_allele_col_name=eqtl_col_dict['ref'],
                                  gbeta_col_name=gwas_col_dict['beta'])
        # logging.info(f'======> Generating LD for significant SNP {range_lead_snp}: {datetime.now()}')
        # # Generate ld file from subset vcf file, plink will append .ld to the --out parameter
        # output_ld_path = os.path.join(coloc_input_dir, f'{range_lead_snp}_chr{chrom}')
        # output_ld_full_path = f'{output_ld_path}.ld'
        # ld_snp_list_path = f'{output_ld_path}.snplist'
        # logging.info(f'ld full path: {output_ld_full_path}')
        # utils.delete_file_if_exists(output_ld_full_path)
        # utils.delete_file_if_exists(ld_snp_list_path)
        # plink_path = gdp.global_config['tool_path']['plink']
        # os.system(f'{plink_path} -vcf {output_vcf_path} --r --matrix  --write-snplist --out {output_ld_path}')
        # if not os.path.exists(output_ld_full_path) or os.path.getsize(output_ld_full_path) <= 0:
        #     logging.info(f'No generated ld file for SNP {range_lead_snp}')
        #     continue

        # Now gwas/eqtl/vcf have the same num of rows, write candidate data to file
        # Reverse GWAS&eQTL column mapping key-value and pass to R so that dataframe in R has fixed column names
        if ('varbeta' not in eqtl_col_dict.keys() or eqtl_col_dict.get('varbeta') is None) and (
                'se' in eqtl_col_dict.keys() and eqtl_col_dict.get('se') is not None):
            eqtl_trait_df['varbeta'] = eqtl_trait_df[eqtl_col_dict['se']] ** 2
        eqtl_trait_df.rename({v: k for k, v in eqtl_col_dict.items()}, axis='columns', inplace=True)
        if ('varbeta' not in gwas_col_dict.keys() or gwas_col_dict.get('varbeta') is None) and (
                'se' in gwas_col_dict.keys() and gwas_col_dict.get('se') is not None):
            candidate_gwas_df['varbeta'] = candidate_gwas_df[gwas_col_dict['se']] ** 2
        candidate_gwas_df.rename({v: k for k, v in gwas_col_dict.items()}, axis='columns', inplace=True)
        eqtl_trait_df.to_csv(coloc_eqtl_input_path, sep=const.output_spliter, header=True, index=False)
        del eqtl_trait_df
        candidate_gwas_df.to_csv(coloc_gwas_input_path, sep=const.output_spliter, header=True, index=False)
        del candidate_gwas_df

        logging.debug(f'Running coloc for significant SNP {range_lead_snp} on ',
                      f'GWAS input file {coloc_gwas_input_path} and eQTL input file {coloc_eqtl_input_path}')
        rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'coloc.R')
        os.system(f'Rscript --no-save --no-restore {rscript_path} '
                  f'{output_file} {coloc_gwas_input_path} {coloc_eqtl_input_path} '
                  f'{gwas_sample_size} {eqtl_sample_size} {gwas_type} {eqtl_type} {p1} {p2} {p12}')
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
        single_result_list = []
        for single_result in os.listdir(output_dir):
            if not (single_result.endswith('.tsv') or single_result.endswith('.tsv.gz')):
                continue
            single_result_list.append(pd.read_table(os.path.join(output_dir, single_result)))
        if len(single_result_list) == 0:
            return
        report_df = pd.concat(single_result_list)
        report_df.sort_values(by=['overall_H4', 'SNP.PP.H4'], ascending=False, inplace=True)
        report_df.drop_duplicates(subset=['snp', 'SNP.PP.H4', 'gene_id'], inplace=True)
        report_df = report_df.round(4)
        report_df.to_csv(final_result_file, sep=const.output_spliter, header=True, index=False)
        ## FDR threshold
        fdrthreshold_outfile = os.path.join(working_dir, 'analyzed', 'fdr_threshold.txt')
        prob_thresh = prob_fdr.calc_threshold_for_prob_rpt(final_result_file, 'overall_H4')
        print(f"colocthreshold: {prob_thresh}")
        with open(fdrthreshold_outfile, 'w') as f:
            f.write(str(prob_thresh))
        f.close()

    def get_output_file(self, working_dir):
        _output_file_name = f'{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        return os.path.join(working_dir, 'analyzed', _output_file_name)


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    if len(os.listdir(glob_processor.gwas_cluster_output_dir)) == 0 or not os.path.exists(
            glob_processor.eqtl_output_report):
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, Coloc.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _eqtl_sample_size = glob_processor.global_config['input']['eqtl'].get('sample_size', 948)
    _gwas_type = glob_processor.global_config['input']['gwas'].get('type', 'cc')
    _eqtl_type = glob_processor.global_config['input']['eqtl'].get('type', 'quant')
    coloc = Coloc()
    coloc.run(_working_dir,
              gdp.Processor.VAR_ID_COL_NAME,
              glob_processor.gwas_cluster_output_dir,
              glob_processor.gwas_col_dict,
              _gwas_sample_size,
              glob_processor.eqtl_output_report,
              glob_processor.eqtl_output_dir,
              glob_processor.eqtl_col_dict,
              _eqtl_sample_size,
              _gwas_type,
              _eqtl_type)
