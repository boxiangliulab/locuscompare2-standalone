import ast
import logging
import os
import re
import sys
from datetime import datetime
from pathlib import Path
import yaml
import pandas as pd
import pyranges as pr

import smr.qtl_data_processor as edp
from fdr import pval_fdr

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import utils, global_data_process as gdp, constants as const


class Smr:
    COLOC_TOOL_NAME = 'smr'

    def __init__(self):
        logging.info('init Smr')

    def run(self,
            working_dir=None,
            var_id_col_name=None,
            genecode_file=None,
            gwas_filter_file=None,
            gwas_chrom_group_dir=None,
            gwas_col_dict=None,
            gwas_sample_size=None,
            qtl_output_report=None,
            qtl_grouped_dir=None,
            qtl_col_dict=None,
            qtl_p_thresh=None,
            subset_vcf_dir=None,
            subset_vcf_file_pattern=None,
            ld_ref_dir=None,
            ld_ref_file_pattern=None,
            tools_config_file=None,
            min_matching_number=0):
        
        start_time = datetime.now()
        logging.info(f'run_smr start at: {start_time}')
        input_dir = self.__get_input_dir(working_dir)
        output_dir = self.__get_output_dir(working_dir)
        Path(input_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        smr_custom_params = utils.get_tools_params(self.COLOC_TOOL_NAME, tools_config_file=tools_config_file,
                                                   params_without_value=['heidi-off'])

        pval_filtered_gwas_df = pd.read_table(gwas_filter_file, sep=const.column_spliter,
                                              usecols=[gwas_col_dict['chrom'], gwas_col_dict['position']],
                                              dtype={gwas_col_dict['chrom']: 'category',
                                                     gwas_col_dict['position']: 'Int64'})
        genecode_df = self.__prepare_genecode(genecode_file)
        qtl_summary_df = pd.read_table(qtl_output_report, sep=const.column_spliter,
                                        dtype={qtl_col_dict['chrom']: 'category'})
        # Loop to process all eQTL trait file
        for chrom_group_file in os.listdir(gwas_chrom_group_dir):
            if not (chrom_group_file.endswith('.tsv') or chrom_group_file.endswith('.tsv.gz')):
                continue
            chr_nums = re.findall(r'\d+', chrom_group_file)
            if len(chr_nums) == 0:
                continue
            chrom = str(chr_nums[0])
            gwas_file = os.path.join(gwas_chrom_group_dir, chrom_group_file)
            gwas_chrom_df = pd.read_table(gwas_file, sep=const.column_spliter,
                                          usecols=[
                                              gwas_col_dict['chrom'],
                                              gwas_col_dict['position'],
                                              gwas_col_dict['effect_allele'],
                                              gwas_col_dict['other_allele'],
                                              gwas_col_dict['beta'],
                                              gwas_col_dict['se'],
                                              gwas_col_dict['pvalue'],
                                              var_id_col_name],
                                          dtype={
                                              gwas_col_dict['chrom']: 'category',
                                              gwas_col_dict['position']: 'Int64',
                                            #   gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                            #   gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)
                                              gwas_col_dict['effect_allele']: 'category',
                                              gwas_col_dict['other_allele']: 'category',
                                            })
            if gwas_chrom_df.shape[0] < min_matching_number:
                continue
            for _, row in qtl_summary_df.iterrows():
                if chrom != str(row.loc['chrom']):
                    continue
                qtl_gene_file = os.path.join(qtl_grouped_dir, chrom, row.loc['gene'], row.loc['pheno_file'])
                pheno_id = utils.get_pheno_name(qtl_gene_file)
                gene_id = row.loc['gene']
                qtl_positions = ast.literal_eval(row.loc['positions'])
                if len(qtl_positions) < min_matching_number:
                    continue
                if not pval_filtered_gwas_df[pval_filtered_gwas_df[gwas_col_dict['chrom']] == chrom].loc[
                       :, gwas_col_dict['position']].isin(qtl_positions).any():
                    logging.debug(
                        f'No intersection between gwas risk SNPs and qtl risk SNPs for phenotype {pheno_id}: {datetime.now()}')
                    continue
                # qtl_trait_df = pd.read_table(qtl_gene_file, sep=const.column_spliter, dtype=qtl_type_dict)
                # 这里过滤不过滤都一样, smr代码会自行进行过滤
                # candidate_qtl_trait_df = utils.filter_data_frame_by_p_value(qtl_trait_df,
                #                                                              global_config['p-value_threshold']['qtl'],
                #                                                              qtl_col_dict['pvalue'], inplace=False)
                # candidate_qtl_trait_df = qtl_trait_df.copy()
                # del qtl_trait_df
                candidate_qtl_trait_df = pd.read_table(qtl_gene_file, sep=const.column_spliter,
                                                        usecols=[
                                                            qtl_col_dict['chrom'],
                                                            qtl_col_dict['position'],
                                                            qtl_col_dict['alt'],
                                                            qtl_col_dict['ref'],
                                                            qtl_col_dict['beta'],
                                                            qtl_col_dict['se'],
                                                            qtl_col_dict['pvalue'],
                                                            var_id_col_name],
                                                        dtype={
                                                            qtl_col_dict['chrom']: 'category',
                                                            qtl_col_dict['position']: 'Int64',
                                                            # qtl_col_dict['alt']: ,
                                                            # qtl_col_dict['ref']: ,
                                                        })
                candidate_gwas_df = gwas_chrom_df.drop(index=gwas_chrom_df[
                    ~gwas_chrom_df[var_id_col_name].isin(candidate_qtl_trait_df[var_id_col_name])].index,
                                       inplace=False)
                if len(candidate_gwas_df) < min_matching_number:
                    print(f"len(candidate_gwas_df) < min_matching_number")
                    continue
                candidate_gwas_df.reset_index(drop=True, inplace=True)
                utils.drop_non_intersect_rows(candidate_qtl_trait_df, var_id_col_name,
                                              candidate_gwas_df, var_id_col_name)
                if len(candidate_gwas_df) < min_matching_number:
                    print(f"len(candidate_gwas_df) < min_matching_number2")
                    continue
                candidate_qtl_trait_df.reset_index(drop=True, inplace=True)
                gwas_input_path = os.path.join(input_dir, f'gwas_chr{chrom}_{gene_id}_{pheno_id}.tsv')
                qtl_input_path = os.path.join(input_dir, f'qtl_chr{chrom}_{gene_id}_{pheno_id}.tsv')
                matching_rpt_name = f'{utils.get_file_name(subset_vcf_file_pattern).format(chrom, gene_id, pheno_id)}.tsv'
                print(f"matching_rpt_name: {matching_rpt_name}")
                vcf_matching_file = os.path.join(subset_vcf_dir, 'matching', matching_rpt_name)
                print(f"vcf_matching_file: {vcf_matching_file}")
                if not os.path.exists(vcf_matching_file) or os.path.getsize(vcf_matching_file) <= 0:
                    logging.warning(f'No generated vcf file for gene {gene_id} {pheno_id}')
                    continue
                vcf_matching_df = pd.read_table(vcf_matching_file,
                                                sep=const.column_spliter, header=0,
                                                usecols=[qtl_col_dict['position'], var_id_col_name],
                                                dtype={qtl_col_dict['position']: 'Int64'})
                if len(vcf_matching_df) < min_matching_number:
                    continue
                # Drop GWAS rows that does not have vcf records
                # SNP in vcf_matching_df is subset of SNP in candidate_gwas_df, so it's fine to drop intersect rows here
                utils.drop_non_intersect_rows(candidate_gwas_df, var_id_col_name, vcf_matching_df, var_id_col_name)
                # Drop eQTL rows that does not have vcf records
                utils.drop_non_intersect_rows(candidate_qtl_trait_df, var_id_col_name, vcf_matching_df,
                                              var_id_col_name)
                del vcf_matching_df
                if len(candidate_gwas_df) < min_matching_number:
                    continue
                candidate_gwas_df.reset_index(drop=True, inplace=True)
                candidate_qtl_trait_df.reset_index(drop=True, inplace=True)
                qtl_flist_file = os.path.join(input_dir, f'{gene_id}_{pheno_id}_flist.txt')
                qtl_besd_file = os.path.join(input_dir, f'{gene_id}_{pheno_id}')
                utils.delete_file_if_exists(qtl_flist_file)
                utils.delete_file_if_exists(qtl_besd_file)
                output_ld_ref_path = os.path.join(ld_ref_dir, ld_ref_file_pattern.format(chrom, gene_id, pheno_id))
                print(f"output_ld_ref_path: {output_ld_ref_path}")
                if not os.path.exists(f'{output_ld_ref_path}.bim') or not os.path.exists(
                        f'{output_ld_ref_path}.bed') or not os.path.exists(f'{output_ld_ref_path}.fam'):
                    logging.warning(f'No generated ld ref file for gene {gene_id} {pheno_id}')
                    continue
                # Adjust allele order and then set GWAS eaf to NA
                # utils.adjust_allele_order(candidate_gwas_df,
                #                           gwas_col_dict['effect_allele'],
                #                           gwas_col_dict['other_allele'],
                #                           gwas_col_dict['chrom'],
                #                           gwas_col_dict['position'],
                #                           candidate_qtl_trait_df,
                #                           ref_df_chrom_col_name=qtl_col_dict['chrom'],
                #                           ref_df_pos_col_name=qtl_col_dict['position'],
                #                           ref_df_alt_allele_col_name=qtl_col_dict['alt'],
                #                           ref_df_ref_allele_col_name=qtl_col_dict['ref'],
                #                           gbeta_col_name=gwas_col_dict['beta'])
                # Keep only SMR recognized COJO format columns
                candidate_gwas_df['eaf'] = 'NA'
                candidate_gwas_df = candidate_gwas_df[
                    [var_id_col_name, gwas_col_dict['effect_allele'],
                     gwas_col_dict['other_allele'],
                     'eaf', gwas_col_dict['beta'], gwas_col_dict['se'],
                     gwas_col_dict['pvalue']]]
                candidate_gwas_df['N'] = gwas_sample_size
                gwas_col = {v: k for k, v in gwas_col_dict.items()}
                # gwas file first col name has to be 'snp'
                gwas_col.update({var_id_col_name: 'snp'})
                candidate_gwas_df.rename(gwas_col, axis='columns', inplace=True)
                candidate_gwas_df.to_csv(gwas_input_path, sep=const.output_spliter, header=True, index=False)
                del candidate_gwas_df
                # eQTL data does not have eaf column, set as NA
                candidate_qtl_trait_df['eaf'] = 'NA'
                candidate_qtl_trait_df = candidate_qtl_trait_df[
                    [qtl_col_dict['chrom'], var_id_col_name, qtl_col_dict['position'],
                     qtl_col_dict['alt'],
                     qtl_col_dict['ref'], 'eaf', qtl_col_dict['beta'], qtl_col_dict['se'],
                     qtl_col_dict['pvalue']]]
                qtl_col = {v: k for k, v in qtl_col_dict.items()}
                # esd file first column name has to be 'Chr', rest col names can be arbitrary strings
                qtl_col.update({qtl_col_dict['chrom']: 'Chr'})
                candidate_qtl_trait_df.rename(qtl_col, axis='columns', inplace=True)
                candidate_qtl_trait_df.to_csv(qtl_input_path, sep=const.output_spliter, header=True, index=False)
                print(f'qtl_input_path: {qtl_input_path}')
                del candidate_qtl_trait_df
                genecode_df_idx = genecode_df[genecode_df['gene_id'] == gene_id].index
                if len(genecode_df_idx) == 0:
                    print(f"len(genecode_df_idx) == 0")
                    continue
                probe_bp = genecode_df.at[genecode_df_idx[0], 'position']
                gene_name = genecode_df.at[genecode_df_idx[0], 'gene_name']
                gene_strand = genecode_df.at[genecode_df_idx[0], 'strand']
                print(f"probe_bp: {probe_bp}")
                print(f"gene_name: {gene_name}")
                print(f"gene_strand: {gene_strand}")

                print(f'qtl_flist_file: {qtl_flist_file}')
                with open(qtl_flist_file, mode='w') as flist_ptr:
                    # flist file first column name has to be 'Chr', can have some extra cols
                    flist_ptr.write(f'Chr\tProbeID\tGeneticDistance\tProbeBp\tGene\tOrientation\tPathOfEsd\n')
                    flist_ptr.write(
                        f'{chrom}\t{gene_id}_{pheno_id}\t0\t{probe_bp}\t{gene_name}\t{gene_strand}\t{qtl_input_path}\n')
                logging.info(f'Generating besd for gene {gene_id}')
                print(f'smr --eqtl-flist {qtl_flist_file} --make-besd --out {qtl_besd_file}')
                os.system(f'smr --eqtl-flist {qtl_flist_file} --make-besd --out {qtl_besd_file}')
                if not os.path.exists(f'{qtl_besd_file}.besd') or not os.path.exists(
                        f'{qtl_besd_file}.esi') or not os.path.exists(f'{qtl_besd_file}.epi'):
                    logging.warning(f'No generated besd file for gene {gene_id}')
                    print(f'No generated besd file for gene {gene_id}')
                    continue
                logging.debug(f'Running SMR on ',
                              f'GWAS input file {gwas_input_path} and eQTL input file {qtl_besd_file}',
                              f' and LD ref file: {output_ld_ref_path}')
                gene_out_result = os.path.join(output_dir, f'{gene_id}_{pheno_id}_out')
                print(f"gene_out_result: {gene_out_result}")
                smr_cmd = f'smr --bfile {output_ld_ref_path} ' \
                          f'--gwas-summary {gwas_input_path} ' \
                          f'--beqtl-summary {qtl_besd_file} ' \
                          f'--peqtl-smr {qtl_p_thresh} ' \
                          f'--out {gene_out_result}'
                cmd_params = ''
                if smr_custom_params and smr_custom_params != '':
                    cmd_params += smr_custom_params
                os.system(f'{smr_cmd} {cmd_params}')
                gene_out_result_file = f'{gene_out_result}.smr'
                if not utils.file_exists(gene_out_result_file):
                    logging.warning(f'SMR: {gene_out_result_file} does not exist')
                    continue
                gene_result_df = pd.read_table(gene_out_result_file)
                gene_result_df.rename(columns={'probeID': 'phenotype_id', 'ProbeChr': 'chrom'}, inplace=True)
                gene_result_df.to_csv(gene_out_result_file, sep=const.output_spliter, header=True, index=False)
        output_file = self.get_output_file(working_dir)
        Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
        self.__analyze_result(output_dir, output_file)

        # fdrthreshold_outfile = os.path.join(working_dir, 'analyzed', 'fdr_threshold.txt')
        if not os.path.exists(output_file) or os.path.getsize(output_file) <= 0:
            ## FDR threshold
            # config = {
            #     'value': 0,
            #     'note': "No result found",
            # }
            # with open(fdrthreshold_outfile, 'w') as file:
            #     yaml.dump(config, file, default_flow_style=False, sort_keys=False)
            logging.warning(f'Process completed, duration {datetime.now() - start_time}, no result found')
        else:
            logging.info(
                f'Process completed, duration {datetime.now() - start_time}, check {output_file} for result!')
            ## FDR threshold
            # pval_thresh, notes = pval_fdr.calc_threshold_for_pval_rpt(output_file, 'p_SMR', working_dir)
            # config = {
            #     'value': float(pval_thresh),
            #     'note': notes,
            # }
            # with open(fdrthreshold_outfile, 'w') as file:
            #     yaml.dump(config, file, default_flow_style=False, sort_keys=False)
        return output_file
    

    def get_output_file(self, working_dir):
        _output_file_name = f'{Smr.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        return os.path.join(working_dir, 'analyzed', _output_file_name)

    def __get_input_dir(self, working_dir):
        return os.path.join(working_dir, 'input')

    def __get_output_dir(self, working_dir):
        return os.path.join(working_dir, 'output')

    def __prepare_genecode(self, genecode_file):
        genecode_df = pr.read_gtf(genecode_file, as_df=True)
        genecode_df.drop(labels=genecode_df[genecode_df['Feature'] != 'gene'].index, inplace=True)
        genecode_df.reset_index(drop=True, inplace=True)
        genecode_df.drop(columns=[col for col in genecode_df if
                                  col not in ['Start', 'End', 'Strand', 'gene_id', 'gene_name']],
                         inplace=True)
        genecode_df.rename({'Start': 'start', 'End': 'end', 'Strand': 'strand'},
                           axis='columns', inplace=True)
        # gene_id_df = genecode_df['gene_id'].str.split('.', n=1, expand=True)
        # genecode_df.loc[:, 'gene_id'] = gene_id_df[0]
        genecode_df.loc[:, 'position'] = pd.NA
        genecode_df['position'].mask(genecode_df['strand'] == '+', genecode_df['start'], inplace=True)
        genecode_df['position'].mask(genecode_df['position'].isna(), genecode_df['end'], inplace=True)
        genecode_df.drop(columns=['start', 'end'], inplace=True)
        return genecode_df

    def __analyze_result(self, output_dir,
                         final_report_file):
        # Merge every single result
        single_result_list = []
        for single_result in os.listdir(output_dir):
            if not single_result.endswith('.smr'):
                continue
            single_result_list.append(pd.read_table(os.path.join(output_dir, single_result)))
        if len(single_result_list) == 0:
            return
        report_df = pd.concat(single_result_list)
        report_df.sort_values(by='p_SMR', ascending=True, inplace=True)
        report_df.drop_duplicates(subset=['topSNP', 'p_SMR', 'phenotype_id'], inplace=True)
        report_df['gene_id'] = report_df['phenotype_id'].apply(lambda x: x.split('_')[0])
        report_df['phenotype_id'] = report_df['phenotype_id'].apply(lambda x: '_'.join(x.split('_')[1:]))
        if len(report_df) <= 0:
            utils.delete_file_if_exists(final_report_file)
        else:
            report_df['p_SMR'] = report_df['p_SMR'].apply(lambda x: "{:.3e}".format(x))
            report_df.to_csv(final_report_file, sep=const.output_spliter, header=True, index=False)


if __name__ == '__main__':
    if len(sys.argv) < 6:
        raise ValueError(f'These arguments are required: \n'
                         f'gwas group file dir\n'
                         f'subset vcf dir\n'
                         f'subset vcf pattern\n'
                         f'ld ref file dir\n'
                         f'ld ref file pattern\n')
    glob_processor = gdp.Processor()
    if not os.path.exists(glob_processor.gwas_filter_file) or os.path.getsize(
            glob_processor.gwas_filter_file) <= 0 or not os.path.exists(
        glob_processor.qtl_output_report) or os.path.getsize(glob_processor.qtl_output_report) <= 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, Smr.COLOC_TOOL_NAME)
    _gwas_chrom_group_dir = sys.argv[1]
    if len(os.listdir(_gwas_chrom_group_dir)) == 0:
        raise ValueError(f'Dependant files not found, did you run gwas_data_processor?')
    _subset_vcf_dir = sys.argv[2]
    _ld_ref_dir = sys.argv[4]
    if len(os.listdir(_ld_ref_dir)) == 0 or len(os.listdir(_subset_vcf_dir)) == 0:
        raise ValueError(f'Dependant files not found, did you run qtl_data_processor?')

    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _qtl_sample_size = glob_processor.global_config['input']['qtl'].get('sample_size', 948)
    _genecode_file = glob_processor.global_config['input']['genecode']
    smr = Smr()
    smr.run(_working_dir,
            gdp.Processor.VAR_ID_COL_NAME,
            _genecode_file,
            glob_processor.gwas_filter_file,
            _gwas_chrom_group_dir,
            glob_processor.gwas_col_dict,
            _gwas_sample_size,
            glob_processor.qtl_output_report,
            glob_processor.qtl_grouped_dir,
            glob_processor.qtl_col_dict,
            glob_processor.config_holder.qtl_p_threshold,
            _subset_vcf_dir,
            sys.argv[3],
            _ld_ref_dir,
            sys.argv[5])
