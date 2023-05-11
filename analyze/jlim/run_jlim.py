# sys.path.append('/Volumes/HD/gitrepo/colocalization-tools')
import ast
import concurrent
import logging
import os
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import pandas as pd

from common import coloc_utils as utils, global_data_process as gdp, constants as const


class Jlim:
    COLOC_TOOL_NAME = 'jlim'

    def __init__(self):
        logging.info('init Jlim')
        # defult use preprocessed eqtl
        self.is_ref_db_eqtlcatalogue = False

        # estimated JLIM p-value
        self.jlim_output_file_pvalue_clo_name = 'pvalue'
        self.shell_command_jlim_execute = \
            '{} --maintr-file {} --maintr-colname-chr {} --maintr-colname-bp {} --maintr-colname-p {}' \
            ' --sectr-file {}{} --ref-ld {} --index-snp {}:{} --output-file {}' \
            ' --min-SNPs-count 5 --sectr-sample-size {} {}'
        self.shell_command_jlim_rscript_run = 'Rscript -e "jlimR:::run.jlim()" ARGSTART $*'
        self.shell_command_jlim_ref_db_eqtlcatalogue = '--sectr-ref-db eQTLCatalogue'
        self.shell_command_jlim_sectr_file_columns = \
            ' --sectr-colname-chr {} --sectr-colname-gene {} --sectr-colname-bp {} --sectr-colname-p {}'

    def run(self,
            gwas_col_dict=None,
            eqtl_col_dict=None,
            var_id_col_name=None,
            eqtl_file_path=None,
            reference_genotype_panel_ld_ref_file=None,
            eqtl_sample_size=None,
            working_dir=None,
            gwas_cluster_output_dir=None,
            eqtl_output_report=None,
            gwas_filter_file=None,
            eqtl_output_dir=None,
            parallel=False):

        start_time = datetime.now()
        logging.info(f'jlim start at: {start_time}')

        # find GWAS/eQTL SNP column name from input file
        gwas_rsid_col_name = var_id_col_name
        eqtl_rsid_col_name = var_id_col_name

        gwas_type_dict = {gwas_col_dict['chrom']: 'category',
                          gwas_col_dict['position']: 'Int64'}

        eqtl_type_dict = {eqtl_col_dict['chrom']: 'category',
                          eqtl_col_dict['position']: 'Int64'}

        output_base_dir = working_dir
        output_gwas_dir = f'{output_base_dir}/gwas'
        output_eqtl_dir = f'{output_base_dir}/eqtl'
        output_jlim_output_dir = f'{output_base_dir}/jlim_output'
        output_analyze_output_dir = f'{output_base_dir}/analyzed'
        output_jlim_output_whole_file = f'{output_analyze_output_dir}/{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}'

        Path(output_jlim_output_dir).mkdir(parents=True, exist_ok=True)

        if self.is_ref_db_eqtlcatalogue:
            self.__run_jlim_ref_db_eQTLCatalogue(gwas_cluster_output_dir=gwas_cluster_output_dir,
                                                 output_jlim_output_dir=output_jlim_output_dir,
                                                 gwas_col_dict=gwas_col_dict,
                                                 eqtl_file_path=eqtl_file_path,
                                                 reference_genotype_panel_ld_ref_file=reference_genotype_panel_ld_ref_file,
                                                 eqtl_sample_size=eqtl_sample_size)
        else:
            self.__run_jlim_preprocessed_eqtl(output_gwas_dir=output_gwas_dir,
                                              output_eqtl_dir=output_eqtl_dir,
                                              eqtl_output_report=eqtl_output_report,
                                              eqtl_col_dict=eqtl_col_dict,
                                              gwas_col_dict=gwas_col_dict,
                                              gwas_type_dict=gwas_type_dict,
                                              eqtl_type_dict=eqtl_type_dict,
                                              gwas_cluster_output_dir=gwas_cluster_output_dir,
                                              gwas_filter_file=gwas_filter_file,
                                              eqtl_output_dir=eqtl_output_dir,
                                              gwas_rsid_col_name=gwas_rsid_col_name,
                                              eqtl_rsid_col_name=eqtl_rsid_col_name,
                                              reference_genotype_panel_ld_ref_file=reference_genotype_panel_ld_ref_file,
                                              eqtl_sample_size=eqtl_sample_size,
                                              output_jlim_output_dir=output_jlim_output_dir,
                                              parallel=parallel)
        reppert_file = self.analyze_result(output_jlim_output_dir, output_analyze_output_dir,
                                           output_jlim_output_whole_file)
        logging.info(f'jlim complete at: {datetime.now()},duration: {datetime.now() - start_time}')
        return reppert_file

    def __run_jlim_ref_db_eQTLCatalogue(self,
                                        gwas_cluster_output_dir=None,
                                        output_jlim_output_dir=None,
                                        gwas_col_dict=None,
                                        eqtl_file_path=None,
                                        reference_genotype_panel_ld_ref_file=None,
                                        eqtl_sample_size=None):

        for gwas_cluster_file_name in os.listdir(gwas_cluster_output_dir):
            if gwas_cluster_file_name.startswith('chr'):
                file_name = gwas_cluster_file_name.split('.')[0]
                jlim_input_file_name_path = f'{gwas_cluster_output_dir}/{gwas_cluster_file_name}'
                output_jlim_file = f'{output_jlim_output_dir}/output_{gwas_cluster_file_name}.out'
                file_name_list = file_name.split('-')[0].split('_')
                index_chr = file_name_list[0].strip("chr")
                index_position = file_name_list[1]
                os.system(self.shell_command_jlim_execute.format(self.shell_command_jlim_rscript_run,
                                                                 jlim_input_file_name_path,
                                                                 gwas_col_dict['chrom'],
                                                                 gwas_col_dict['position'],
                                                                 gwas_col_dict['pvalue'],
                                                                 eqtl_file_path, '',
                                                                 reference_genotype_panel_ld_ref_file, index_chr,
                                                                 index_position,
                                                                 output_jlim_file,
                                                                 eqtl_sample_size,
                                                                 self.shell_command_jlim_ref_db_eqtlcatalogue))

    def __run_jlim_preprocessed_eqtl(self,
                                     output_gwas_dir=None,
                                     output_eqtl_dir=None,
                                     eqtl_output_report=None,
                                     eqtl_col_dict=None,
                                     gwas_col_dict=None,
                                     gwas_type_dict=None,
                                     eqtl_type_dict=None,
                                     gwas_cluster_output_dir=None,
                                     gwas_filter_file=None,
                                     eqtl_output_dir=None,
                                     gwas_rsid_col_name=None,
                                     eqtl_rsid_col_name=None,
                                     output_jlim_output_dir=None,
                                     reference_genotype_panel_ld_ref_file=None,
                                     eqtl_sample_size=None,
                                     parallel=False):
        start_time = datetime.now()

        Path(output_gwas_dir).mkdir(parents=True, exist_ok=True)
        Path(output_eqtl_dir).mkdir(parents=True, exist_ok=True)

        eqtl_summary_df = pd.read_csv(eqtl_output_report, sep=const.column_spliter,
                                      dtype={eqtl_col_dict['chrom']: 'category'})
        pval_filtered_gwas_df = pd.read_table(gwas_filter_file, sep=const.column_spliter,
                                              dtype=gwas_type_dict)

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

        gwas_chroms = pval_filtered_gwas_df[gwas_col_dict['chrom']].unique().tolist()

        if parallel:
            with ThreadPoolExecutor(max_workers=10) as executor:
                futures = []
                for _, row in eqtl_summary_df.iterrows():
                    chrom = str(row.loc['chrom'])
                    if chrom not in gwas_chroms:
                        return
                    eqtl_gene_file = os.path.join(eqtl_output_dir, chrom, row.loc['gene_file'])
                    gene_id = utils.get_file_name(eqtl_gene_file)
                    for gwas_range_file in gwas_range_files[chrom]:
                        futures.append(executor.submit(self.process_gene,
                                                       gwas_range_file, gwas_type_dict, gwas_col_dict, row,
                                                       eqtl_gene_file, eqtl_type_dict,
                                                       eqtl_rsid_col_name, gwas_rsid_col_name, output_gwas_dir,
                                                       output_eqtl_dir,
                                                       gene_id, eqtl_col_dict,
                                                       output_jlim_output_dir, chrom,
                                                       reference_genotype_panel_ld_ref_file,
                                                       eqtl_sample_size))

                for future in concurrent.futures.as_completed(futures):
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error('Get result generated an exception: %s' % exc)
        else:
            for _, row in eqtl_summary_df.iterrows():
                chrom = str(row.loc['chrom'])
                if chrom not in gwas_chroms:
                    return
                eqtl_gene_file = os.path.join(eqtl_output_dir, chrom, row.loc['gene_file'])
                gene_id = utils.get_file_name(eqtl_gene_file)
                for gwas_range_file in gwas_range_files[chrom]:
                    self.process_gene(gwas_range_file, gwas_type_dict, gwas_col_dict, row,
                                      eqtl_gene_file, eqtl_type_dict,
                                      eqtl_rsid_col_name, gwas_rsid_col_name, output_gwas_dir,
                                      output_eqtl_dir,
                                      gene_id, eqtl_col_dict,
                                      output_jlim_output_dir, chrom,
                                      reference_genotype_panel_ld_ref_file,
                                      eqtl_sample_size)

        logging.info(f'Process completed, duration {datetime.now() - start_time} seconds')

    def process_gene(self, gwas_range_file, gwas_type_dict, gwas_col_dict, row, eqtl_gene_file, eqtl_type_dict,
                     eqtl_rsid_col_name, gwas_rsid_col_name, output_gwas_dir, output_eqtl_dir, gene_id, eqtl_col_dict,
                     output_jlim_output_dir, chrom, reference_genotype_panel_ld_ref_file, eqtl_sample_size):
        range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
        candidate_gwas_df = pd.read_table(gwas_range_file, sep=const.column_spliter, dtype=gwas_type_dict)
        if len(candidate_gwas_df) <= 1:
            return
        if not candidate_gwas_df.loc[:, gwas_col_dict['position']].isin(
                ast.literal_eval(row.loc['positions'])).any():
            # print(f'======> No intersection between gwas {gwas_range_file} and eqtl {eqtl_gene_file}: {time.time()}')
            return
        eqtl_trait_df = pd.read_table(eqtl_gene_file, sep=const.column_spliter, dtype=eqtl_type_dict)
        candidate_eqtl_trait_df = eqtl_trait_df[
            eqtl_trait_df[eqtl_rsid_col_name].isin(candidate_gwas_df[gwas_rsid_col_name])].copy()
        del eqtl_trait_df
        if len(candidate_eqtl_trait_df) <= 1:
            return
        utils.drop_non_intersect_rows(candidate_eqtl_trait_df, eqtl_rsid_col_name,
                                      candidate_gwas_df, gwas_rsid_col_name)
        if len(candidate_gwas_df) <= 1:
            return
        # if file name is not in file，find the min pvalue of snp to file name
        snp_pos = range_lead_snp.split('_')[1]
        if not candidate_gwas_df[gwas_col_dict['position']].isin([snp_pos]).any():
            target_row = candidate_gwas_df[candidate_gwas_df[gwas_col_dict['position']] == candidate_gwas_df[
                gwas_col_dict['position']].min()].iloc[0]
            logging.debug(f'orign gwas file: {range_lead_snp}')
            range_lead_snp = f'chr{target_row[gwas_col_dict["chrom"]]}_{target_row[gwas_col_dict["position"]]}'
            logging.debug(f'current range_lead_snp:{range_lead_snp}')

        jlim_gwas_input_path = os.path.join(output_gwas_dir, f'gwas_{range_lead_snp}.txt')
        jlim_eqtl_input_path = os.path.join(output_eqtl_dir, f'eqtl_{gene_id}_{range_lead_snp}')

        candidate_eqtl_trait_df['gene_id'] = gene_id
        candidate_eqtl_trait_df.rename({v: k for k, v in eqtl_col_dict.items()}, axis='columns',
                                       inplace=True)
        candidate_gwas_df.rename({v: k for k, v in gwas_col_dict.items()}, axis='columns',
                                 inplace=True)

        candidate_eqtl_trait_df.to_csv(jlim_eqtl_input_path, sep=const.output_spliter, header=True, index=False)
        candidate_gwas_df.to_csv(jlim_gwas_input_path, sep=const.output_spliter, header=True, index=False)
        os.system('gzip -k -f {}'.format(jlim_eqtl_input_path))

        output_jlim_file = f'{output_jlim_output_dir}/output_{gene_id}_{range_lead_snp}.out'
        shell_command_sectr_file_columns_str = self.shell_command_jlim_sectr_file_columns.format(
            'chrom', 'gene_id', 'position', 'pvalue')

        str_com = self.shell_command_jlim_execute.format(self.shell_command_jlim_rscript_run,
                                                         jlim_gwas_input_path,
                                                         'chrom',
                                                         'position',
                                                         'pvalue',
                                                         f'{jlim_eqtl_input_path}.gz',
                                                         shell_command_sectr_file_columns_str,
                                                         reference_genotype_panel_ld_ref_file,
                                                         chrom,
                                                         range_lead_snp.split('_')[1],
                                                         output_jlim_file, eqtl_sample_size,
                                                         '')
        os.system(str_com)

        if utils.file_exists(output_jlim_file):
            # report file add chromesome columns, and drop NAN data
            if os.path.getsize(output_jlim_file):
                output_jlim_df = pd.read_csv(output_jlim_file, sep=const.column_spliter).dropna(axis=0,
                                                                                                how='any')
                if len(output_jlim_df) > 0:
                    output_jlim_df['chrom'] = chrom
                    jlim_output_columns_position = ' actualIdxBP'
                    merge_pd = output_jlim_df
                    be_merge_snp_df = None
                    if gwas_col_dict.get('snp') is not None:
                        be_merge_snp_df = candidate_gwas_df[['position', 'snp']]
                    elif eqtl_col_dict.get('snp') is not None:
                        be_merge_snp_df = candidate_eqtl_trait_df[['position', 'snp']]

                    if be_merge_snp_df is not None:
                        merge_pd = pd.merge(left=output_jlim_df,
                                            right=be_merge_snp_df,
                                            left_on=jlim_output_columns_position,
                                            right_on='position',
                                            how='left').rename(
                            columns={'snp': 'rsid', 'sectrGeneName': 'gene_id'})
                    merge_pd.to_csv(output_jlim_file, header=True, index=False, sep=const.output_spliter)
                else:
                    utils.delete_file_if_exists(output_jlim_file)
            else:
                utils.delete_file_if_exists(output_jlim_file)

        logging.info(f'{output_jlim_file} completed')
        return output_jlim_file

    def analyze_result(self, output_dir, output_analyze_output_dir, output_jlim_output_whole_file):
        jlim_out_list = []
        for jlim_output_file_name in os.listdir(output_dir):
            if jlim_output_file_name.startswith('output'):
                jlim_out_list.append(pd.read_csv(f'{output_dir}/{jlim_output_file_name}', sep=const.column_spliter))

        # Consolidate files and sort pvlaue
        Path(output_analyze_output_dir).mkdir(parents=True, exist_ok=True)
        reppert_file = f'{output_jlim_output_whole_file}.tsv.gz'
        if len(jlim_out_list) > 0:
            pd.concat(jlim_out_list).sort_values(self.jlim_output_file_pvalue_clo_name).to_csv(
                reppert_file, header=True,
                index=False, sep=const.output_spliter)

        return reppert_file


if __name__ == '__main__':
    processor = gdp.Processor()
    global_config = processor.global_config
    eqtl_file_path = f'{global_config["input"]["eqtl"]["file"]}'
    reference_genotype_panel_ld_ref_file = global_config['input']['reference_genotype_panel_ld_ref_file']
    eqtl_sample_size = global_config['input']['eqtl']['sample_size']
    _working_dir = os.path.join(processor.tool_parent_dir, Jlim.COLOC_TOOL_NAME)
    jlim = Jlim()
    jlim.run(gwas_col_dict=processor.gwas_col_dict,
             eqtl_col_dict=processor.eqtl_col_dict,
             var_id_col_name=gdp.Processor.VAR_ID_COL_NAME,
             eqtl_file_path=eqtl_file_path,
             reference_genotype_panel_ld_ref_file=reference_genotype_panel_ld_ref_file,
             eqtl_sample_size=eqtl_sample_size,
             working_dir=_working_dir,
             gwas_cluster_output_dir=processor.gwas_cluster_output_dir,
             eqtl_output_report=processor.eqtl_output_report,
             gwas_filter_file=processor.gwas_filter_file,
             eqtl_output_dir=processor.eqtl_output_dir)
