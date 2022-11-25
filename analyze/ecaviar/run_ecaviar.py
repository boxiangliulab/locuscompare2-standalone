import asyncio
import logging
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import coloc_utils as utils, constants as const


class ECaviar:
    ECAVIAR_TOOL_NAME = 'ecaviar'
    shell_command_run_finemap = 'finemap --sss --in-files {}'

    def __init__(self):
        logging.info('init eCaviar')

    async def run(self, working_dir, candidate_data_dir, gwas_cluster_dir, gwas_preprocessed_file, eqtl_grouped_dir,
                  var_id_col_name,
                  gwas_col_dict, eqtl_col_dict):
        logging.info('eCaviar start to run...')
        Path(f'{working_dir}/analyzed').mkdir(parents=True, exist_ok=True)

        tasks = []
        finemap_snp_files = []
        for gene_dir in os.listdir(candidate_data_dir):
            if not gene_dir.startswith('.'):
                for causal_snp_dir in os.listdir(f'{candidate_data_dir}/{gene_dir}'):
                    if not causal_snp_dir.startswith('.'):
                        candidate_dir = f'{candidate_data_dir}/{gene_dir}/{causal_snp_dir}'
                        finemap_file_name = f'{causal_snp_dir}_{gene_dir}'
                        run_finemap_cmd = self.shell_command_run_finemap.format(
                            f'{candidate_dir}/{finemap_file_name}.in')

                        tasks.append(asyncio.create_task(utils.async_run_cmd(run_finemap_cmd)))
                        finemap_snp_files.append((f'{candidate_dir}/gwas_{finemap_file_name}.snp',
                                                  f'{candidate_dir}/eqtl_{finemap_file_name}.snp'))

        await asyncio.gather(*tasks, return_exceptions=True)

        logging.info(f'finish eCaviar process.')
        logging.info(f'generate eCaviar report.')
        return self.generate_report(working_dir=working_dir, finemap_reports=finemap_snp_files,
                                    gwas_cluster_dir=gwas_cluster_dir,
                                    gwas_preprocessed_file=gwas_preprocessed_file,
                                    eqtl_grouped_dir=eqtl_grouped_dir, var_id_col_name=var_id_col_name,
                                    gwas_col_dict=gwas_col_dict, eqtl_col_dict=eqtl_col_dict)

    def generate_report(self, working_dir, finemap_reports, gwas_cluster_dir, gwas_preprocessed_file, eqtl_grouped_dir,
                        var_id_col_name, gwas_col_dict, eqtl_col_dict):
        output_report_path = f'{working_dir}/analyzed/report.tsv'
        var_ids, chroms, gene_ids, gwas_paths, eqtl_paths, gwas_pips, eqtl_pips, gene_clpps = [], [], [], [], [], [] ,[] ,[]
        for gwas_snp, eqtl_snp in finemap_reports:
            if not utils.file_exists(gwas_snp) or not utils.file_exists(eqtl_snp):
                with open(f'{working_dir}/analyzed/error.log', mode='a') as err_report:
                    err_report.write(f'{gwas_snp} or {eqtl_snp} is not found.\n')
                    continue

            gwas_snp_df = pd.read_csv(gwas_snp, sep=' ', usecols=['snp', 'snp_prob'])
            eqtl_snp_df = pd.read_csv(eqtl_snp, sep=' ', usecols=['snp', 'snp_prob'])
            merge_snp_pd = pd.merge(left=gwas_snp_df, right=eqtl_snp_df, how='inner', on='snp',
                                    suffixes=('_gwas', '_eqtl'))
            merge_snp_pd['snp_clpp'] = merge_snp_pd['snp_prob_gwas'] * merge_snp_pd['snp_prob_eqtl']
            gene_clpp_sum = merge_snp_pd['snp_clpp'].sum()
            gwas_clpp_sum = merge_snp_pd['snp_prob_gwas'].sum()
            eqtl_clpp_sum = merge_snp_pd['snp_prob_eqtl'].sum()
            # gwas snp file name example: gwas_chr11_103802549_ENSG00000170962.snp
            gwas_snp_file = gwas_snp.split('/')[-1]
            var_id = '_'.join(gwas_snp_file.split('_')[1:3])
            chrom = var_id.split('_')[0]
            chrom_num = chrom.replace('chr', '')
            gene_id = gwas_snp_file.split('_')[3].split('.')[0]
            eqtl_path = f'{eqtl_grouped_dir}/{chrom_num}/{gene_id}.tsv'

            var_ids.append(var_id)
            chroms.append(chrom_num)
            gene_ids.append(gene_id)
            gwas_paths.append(gwas_preprocessed_file)
            eqtl_paths.append(eqtl_path)
            gwas_pips.append(gwas_clpp_sum)
            eqtl_pips.append(eqtl_clpp_sum)
            gene_clpps.append(gene_clpp_sum)

        report_df = pd.DataFrame({'var_id': var_ids,
                                  'chrom': chroms,
                                  'gene_id': gene_ids,
                                  'gwas_path': gwas_paths,
                                  'eqtl_path': eqtl_paths,
                                  'gwas_pip': gwas_pips,
                                  'eqtl_pip': eqtl_pips,
                                  'clpp': gene_clpps})

        report_df = utils.mapping_var_id_to_rsid(result_df=report_df,
                                                 result_df_var_id_col_name='var_id',
                                                 result_df_gene_id_col_name='gene_id',
                                                 gwas_preprocessed_file=gwas_preprocessed_file,
                                                 ref_var_id_col_name=var_id_col_name,
                                                 gwas_col_dict=gwas_col_dict,
                                                 eqtl_col_dict=eqtl_col_dict)

        report_df.to_csv(output_report_path, sep=const.column_spliter, index=False)
        return output_report_path


if __name__ == '__main__':
    eCaviar = ECaviar()
    _gwas_col_dict = {'snp': 'oldID', 'chrom': 'hm_chrom', 'position': 'hm_pos', 'beta': 'hm_beta',
                                        'effect_allele': 'hm_effect_allele', 'other_allele': 'hm_other_allele', 'pvalue': 'p_value',
                                        'se': 'standard_error', 'eaf': 'hm_effect_allele_frequency'}
    _eqtl_col_dict = {'snp': 'rsid', 'chrom': 'chromosome', 'position': 'position', 'beta': 'beta', 'alt': 'alt',
                                        'ref': 'ref', 'pvalue': 'pvalue', 'se': 'se', 'gene_id': 'molecular_trait_id', 'maf': 'maf'}
    eCaviar.generate_report(working_dir='/Users/nicklin/Desktop/CAD/report', finemap_reports=None,
                            gwas_cluster_dir='/Volumes/HD/biodata/colocalization-tools/preprocessed/gwas/cad/clustered',
                            gwas_preprocessed_file='/Volumes/HD/biodata/colocalization-tools/preprocessed/gwas/cad/preprocessed.tsv',
                            eqtl_grouped_dir='/Volumes/HD/biodata/colocalization-tools/preprocessed/eqtl/Artery_Aorta/grouped',
                            var_id_col_name='var_id',
                            gwas_col_dict=_gwas_col_dict, eqtl_col_dict=_eqtl_col_dict)
