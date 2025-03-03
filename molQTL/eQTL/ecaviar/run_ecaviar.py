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
        logging.info('init eCaviar')

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
        total_len = len(os.listdir(candidate_data_dir))
        ix = 1
        for gene_dir in os.listdir(candidate_data_dir):
            print(f'gene_dir: {gene_dir}')
            if whether_schedual:
                outputschedule(rownum=ix,
                    totalnum=total_len,
                    current_analysis_order = current_analysis_order,
                    total_numof_analyses=total_numof_analyses,
                    rank_dir=rank_dir)
                ix = ix + 1
            print(f"ix: {ix}")
            if gene_dir.startswith('.'):
                continue

            snp_dirs = os.listdir(os.path.join(candidate_data_dir, gene_dir))
            if not snp_dirs:
                raise ValueError(f"No SNP directories found in {os.path.join(candidate_data_dir, gene_dir)}")

            for causal_snp_dir in os.listdir(os.path.join(candidate_data_dir, gene_dir)):
                print(f'ecaviar processing: {causal_snp_dir}')
                if causal_snp_dir.startswith('.'):
                    print(f"finemap causal_snp_dir: {causal_snp_dir}")
                    continue
                candidate_dir = os.path.join(candidate_data_dir, gene_dir, causal_snp_dir)
                finemap_file_name = f'{causal_snp_dir}_{gene_dir}'
                finemap_file_path = os.path.join(candidate_dir, f'{finemap_file_name}.in')
                if Path(finemap_file_path).exists():
                    run_finemap_cmd = self.shell_command_run_finemap.format(
                        finemap_file_path, finemap_params)
                    logging.info(f"run_finemap_cmd: {run_finemap_cmd}")

                    finemap_snp_files.append((f'{candidate_dir}/gwas_{finemap_file_name}.snp',
                                            f'{candidate_dir}/qtl_{finemap_file_name}.snp'))
                    os.system(run_finemap_cmd)
                else:
                    logging.info(f"{finemap_file_name} not exist")

        logging.info(f'finish eCaviar process.')
        logging.info(f'generate eCaviar report.')
        lst = custom_params.split(' ')  # 将字符串拆分为列表
        if '--n-causal-snps' in lst:
            index = lst.index('--n-causal-snps')  # 获取 '--n' 在列表中的索引
            n_causal_snps = lst[index+1]
        else:
            n_causal_snps = 5
        return self.generate_report(working_dir=working_dir, finemap_reports=finemap_snp_files, n_causal_snps=n_causal_snps)

    def generate_report(self, working_dir, finemap_reports, n_causal_snps):
        output_report_path = f'{working_dir}/analyzed/{self.ECAVIAR_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        var_ids, chroms, phenotype_ids, gwas_pips, eqtl_pips, gene_clpps = [], [], [], [], [], [],
        report_df = pd.DataFrame()
        for gwas_snp, eqtl_snp in finemap_reports:
            if not utils.file_exists(gwas_snp) or not utils.file_exists(eqtl_snp):
                logging.warning(f'{gwas_snp} or {eqtl_snp} is not found.')
                continue
            try:
                # logging.info("generate report for ecaviar")
                gwas_snp_df = pd.read_csv(gwas_snp, sep=' ', usecols=['rsid', 'allele1', 'allele2','prob'])
                eqtl_snp_df = pd.read_csv(eqtl_snp, sep=' ', usecols=['rsid', 'allele1', 'allele2','prob'])
            except Exception as e:
                logging.error(f'Read snp file {gwas_snp} or {eqtl_snp} failed. Error: {e}')
                continue
            
            gwas_snp_df['variant_id'] = 'chr'+gwas_snp_df['rsid']+'_'+gwas_snp_df['allele2']+'_'+gwas_snp_df['allele1']
            eqtl_snp_df['variant_id'] = 'chr'+eqtl_snp_df['rsid']+'_'+eqtl_snp_df['allele2']+'_'+eqtl_snp_df['allele1']

            merge_snp_pd = pd.merge(left=gwas_snp_df, right=eqtl_snp_df, how='inner', on='variant_id',
                                    suffixes=('_gwas', '_eqtl'))
            merge_snp_pd["variant_id"] = merge_snp_pd["variant_id"].str.replace(":", "_")
            merge_snp_pd['snp_clpp'] = merge_snp_pd['prob_gwas'] * merge_snp_pd['prob_eqtl']
            merge_snp_pd = merge_snp_pd.sort_values('snp_clpp', ascending=False)
            merge_snp_pd = merge_snp_pd.iloc[:int(n_causal_snps),:]
            # gene_clpp_sum = merge_snp_pd['snp_clpp'].sum()
            # gwas_clpp_sum = merge_snp_pd['prob_gwas'].sum()
            # eqtl_clpp_sum = merge_snp_pd['prob_eqtl'].sum()
            # gwas_snp file name example: gwas_ENSG00000122026_phenotype_id_chr13_26207391_G_A.snp
            
            gwas_snp_file = gwas_snp.split('/')[-1]
            logging.info(f"gwas_snp_file: {gwas_snp_file}")
            var_id = '_'.join(gwas_snp_file.split('_')[1:5])
            chrom = var_id.split('_')[0]
            chrom_num = chrom.replace('chr', '')
            phenotype_id = '.'.join(gwas_snp_file.split('_')[5].split('.')[:2])
            merge_snp_pd['chrom'] = chrom_num
            merge_snp_pd['lead_variant'] = var_id
            merge_snp_pd['gene_id'] = phenotype_id
            merge_snp_pd = merge_snp_pd[['chrom','lead_variant','variant_id','gene_id','prob_gwas','prob_eqtl','snp_clpp']]
            merge_snp_pd.columns = ['chrom','lead_variant','causal_variant','gene_id','gwas_pip','eqtl_pip','clpp']

            report_df = pd.concat([report_df,merge_snp_pd])
            # var_ids.append(var_id)
            # chroms.append(chrom_num)
            # phenotype_ids.append(phenotype_id)
            # gwas_pips.append(gwas_clpp_sum)
            # eqtl_pips.append(eqtl_clpp_sum)
            # gene_clpps.append(gene_clpp_sum)

        # report_df = pd.DataFrame({'var_id': var_ids,
        #                           'chrom': chroms,
        #                           'gene_id': phenotype_ids,
        #                           'gwas_pip': gwas_pips,
        #                           'eqtl_pip': eqtl_pips,
        #                           'clpp': gene_clpps})
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

        logging.info(f"ecaviar output_report_path: {output_report_path}")
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
