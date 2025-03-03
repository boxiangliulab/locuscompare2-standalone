import logging
import os
import sys
from datetime import datetime
from pathlib import Path
import yaml
import pandas as pd
from fdr import prob_fdr
sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import utils, constants as const

def outputschedule(rownum, totalnum, current_analysis_order, total_numof_analyses, rank_dir):
    calculated_schedule = int(40/total_numof_analyses+ rownum/totalnum * 40/total_numof_analyses + 80/total_numof_analyses * (current_analysis_order - 1))
    logging.warning(f"process ecaviar schedual: {calculated_schedule}")
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()

class ECaviar:
    ECAVIAR_TOOL_NAME = 'ecaviar'
    shell_command_run_finemap = 'finemap --sss --in-files {} {} --log'

    def __init__(self):
        logging.info('init sQTL eCaviar')

    async def run(self, working_dir, candidate_data_dir, parallel=False, 
                  tools_config=None, parallel_worker_num=3, rank_dir = None,
                  current_analysis_order = None, total_numof_analyses = None, whether_schedual = False):

        logging.info('eCaviar start to run...')
        Path(f'{working_dir}/analyzed').mkdir(parents=True, exist_ok=True)
        coros = []
        finemap_snp_files = []
        finemap_params = ''
        custom_params = utils.get_tools_params('ecaviar', tools_config_file=tools_config)
        if custom_params and custom_params != '':
            finemap_params = custom_params
        print(f"finemap_params: {finemap_params}")
        total_len = len(os.listdir(candidate_data_dir))
        ix = 1
        for gene_dir in os.listdir(candidate_data_dir):
            if whether_schedual:
                outputschedule(rownum=ix,
                    totalnum=total_len,
                    current_analysis_order = current_analysis_order,
                    total_numof_analyses=total_numof_analyses,
                    rank_dir=rank_dir)
                ix = ix + 1
            if gene_dir.startswith('.'):
                continue
            for pheno_dir in os.listdir(f'{os.path.join(candidate_data_dir, gene_dir)}'):
                for causal_snp_dir in os.listdir(f'{os.path.join(candidate_data_dir, gene_dir, pheno_dir)}'):
                    if causal_snp_dir.startswith('.'):
                        continue
                    candidate_dir = os.path.join(candidate_data_dir, gene_dir, pheno_dir, causal_snp_dir)
                    print(f"candidate_dir: {candidate_dir}")
                    finemap_file_name = f'{gene_dir}_{pheno_dir}_{causal_snp_dir}'
                    print(f"finemap_file_name: {finemap_file_name}")
                    run_finemap_cmd = self.shell_command_run_finemap.format(
                        f'{candidate_dir}/{finemap_file_name}.in', finemap_params)
                    print(f"ecaviar/finemap cmd: {run_finemap_cmd}")
                    os.system(run_finemap_cmd)
                    # coros.append(utils.async_run_cmd(run_finemap_cmd))
                    # finemap_snp_files.append((f'{candidate_dir}/gwas_{finemap_file_name}.snp',
                    #                         f'{candidate_dir}/qtl_{finemap_file_name}.snp'))

                    # async def run_task():
                    #     try:
                    #         await utils.async_run_cmd(run_finemap_cmd)
                    #         finemap_snp_files.append((f'{candidate_dir}/gwas_{finemap_file_name}.snp',
                    #                                   f'{candidate_dir}/qtl_{finemap_file_name}.snp'))
                    #     except Exception as e:
                    #         logging.error(f"Error running sQTL finemap for {finemap_file_name}: {e}")

                    # Submit each finemap task
                    # coros.append(run_task())  


        # await utils.gather_with_limit(parallel_worker_num if parallel else 1, *coros)
        logging.info(f'finish eCaviar process.')
        logging.info(f'generate eCaviar report.')
        return self.generate_report(working_dir=working_dir, finemap_reports=finemap_snp_files)

    def generate_report(self, working_dir, finemap_reports):
        output_report_path = f'{working_dir}/analyzed/{self.ECAVIAR_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        var_ids, chroms, phenotype_ids, gene_ids, gwas_pips, sqtl_pips, gene_clpps = [], [], [], [], [], [], [],
        for gwas_snp, sqtl_snp in finemap_reports:
            if not utils.file_exists(gwas_snp) or not utils.file_exists(sqtl_snp):
                logging.warning(f'{gwas_snp} or {sqtl_snp} is not found.')
                continue
            try:
                # print("generate report for ecaviar")
                gwas_snp_df = pd.read_csv(gwas_snp, sep=' ', usecols=['rsid', 'allele1', 'allele2','prob'])
                sqtl_snp_df = pd.read_csv(sqtl_snp, sep=' ', usecols=['rsid', 'allele1', 'allele2','prob'])
            except Exception as e:
                logging.error(f'Read snp file {gwas_snp} or {sqtl_snp} failed. Error: {e}')
                continue

            gwas_snp_df['variant_id'] = 'chr'+gwas_snp_df['rsid']+'_'+gwas_snp_df['allele2']+'_'+gwas_snp_df['allele1']
            sqtl_snp_df['variant_id'] = 'chr'+sqtl_snp_df['rsid']+'_'+sqtl_snp_df['allele2']+'_'+sqtl_snp_df['allele1']

            merge_snp_pd = pd.merge(left=gwas_snp_df, right=sqtl_snp_df, how='inner', on='snp',
                                    suffixes=('_gwas', '_sqtl'))
            merge_snp_pd['snp_clpp'] = merge_snp_pd['prob_gwas'] * merge_snp_pd['prob_sqtl']
            gene_clpp_sum = merge_snp_pd['snp_clpp'].sum()
            gwas_clpp_sum = merge_snp_pd['prob_gwas'].sum()
            sqtl_clpp_sum = merge_snp_pd['prob_sqtl'].sum()
            # gwas_snp file name example: ecaviar/candidate/ENSG00000118197.14/chr1_200649226_200659027_clu_2350_-/chr1_201267984_C_A/gwas_ENSG00000118197.14_chr1_200649226_200659027_clu_2350_-_chr1_201267984_C_A.snp
            gwas_snp_file = gwas_snp.split('/')[-1]
            var_id = '_'.join(gwas_snp_file.split('_')[-4:]).split('.')[0]
            chrom = var_id.split('_')[0]
            chrom_num = chrom.replace('chr', '')
            gene_id = gwas_snp_file.split('_')[1]
            phenotype_id = '_'.join(gwas_snp_file.split('_')[2:8])

            var_ids.append(var_id)
            chroms.append(chrom_num)
            gene_ids.append(gene_id)
            phenotype_ids.append(phenotype_id)
            gwas_pips.append(gwas_clpp_sum)
            sqtl_pips.append(sqtl_clpp_sum)
            gene_clpps.append(gene_clpp_sum)

        report_df = pd.DataFrame({'var_id': var_ids,
                                  'chrom': chroms,
                                  'gene_id': gene_ids, 
                                  'phenotype_id': phenotype_ids,
                                  'gwas_pip': gwas_pips,
                                  'sqtl_pip': sqtl_pips,
                                  'clpp': gene_clpps})
        report_df.sort_values(by='clpp', ascending=False, inplace=True)
        # fdrthreshold_outfile = os.path.join(working_dir, 'analyzed', 'fdr_threshold.txt')
        if report_df.shape[0] > 0:
            report_df = report_df.round(4)
            report_df.to_csv(output_report_path, sep=const.column_spliter, index=False)
            # # FDR threshold
            # config = {
            #     'value': 0.01,
            #     'note': "General threshold",
            # }
            # with open(fdrthreshold_outfile, 'w') as file:
            #     yaml.dump(config, file, default_flow_style=False, sort_keys=False)
        # else:
        #     config = {
        #         'value': 1,
        #         'note': "No result found",
        #     }
        #     with open(fdrthreshold_outfile, 'w') as file:
        #         yaml.dump(config, file, default_flow_style=False, sort_keys=False)


        return output_report_path


if __name__ == '__main__':
    eCaviar = ECaviar()
    _gwas_col_dict = {'snp': 'oldID', 'chrom': 'hm_chrom', 'position': 'hm_pos', 'beta': 'hm_beta',
                      'effect_allele': 'hm_effect_allele', 'other_allele': 'hm_other_allele', 'pvalue': 'p_value',
                      'se': 'standard_error'}
    _qtl_col_dict = {'snp': 'rsid', 'chrom': 'chromosome', 'position': 'position', 'beta': 'beta', 'alt': 'alt',
                      'ref': 'ref', 'pvalue': 'pvalue', 'se': 'se', 'phenotype_id': 'molecular_trait_id', 'maf': 'maf'}

    snp_reports = []
    _dir = '/Volumes/HD/nlrp/test_run/processed/default/Artery_Aorta/cad/ecaviar/candidate'
    for sub_dir in os.listdir(_dir):
        if not sub_dir.startswith('.'):
            for snp_dir in os.listdir(f'{_dir}/{sub_dir}'):
                if not snp_dir.startswith('.'):
                    g_rpt = ''
                    e_rpt = ''
                    for snp_rpt in os.listdir(f'{_dir}/{sub_dir}/{snp_dir}'):
                        if snp_rpt.startswith('gwas') and snp_rpt.endswith('.snp'):
                            g_rpt = f'{_dir}/{sub_dir}/{snp_dir}/{snp_rpt}'
                        if snp_rpt.startswith('eqtl') and snp_rpt.endswith('.snp'):
                            e_rpt = f'{_dir}/{sub_dir}/{snp_dir}/{snp_rpt}'
                snp_reports.append((g_rpt, e_rpt))

    eCaviar.generate_report(working_dir='/Volumes/HD/nlrp/test_run', finemap_reports=snp_reports)
