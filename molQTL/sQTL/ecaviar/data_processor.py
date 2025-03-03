import concurrent
import json
import logging
import os.path
import re
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import constants as const, utils, global_data_process as gdp, config

def outputschedule(rownum, totalnum,current_analysis_order, total_numof_analyses, rank_dir):
    calculated_schedule = int(rownum/totalnum * 40/total_numof_analyses + 80/total_numof_analyses * (current_analysis_order - 1))
    print(f"process ecaviar schedual: {calculated_schedule}")
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()

class ECaviarDataProcessor:
    ECAVIAR_TOOL_NAME = 'ecaviar'
    zscore_col = 'zscore'
    shell_command_plink_execute = 'plink --silent --vcf {} --r2 --matrix --mac 1 --write-snplist --out {}'

    def __init__(self):
        logging.info('init sQTL ECaviarDataProcessor')

    def prepare(self, working_dir, gwas_cluster_dir, gwas_cluster_summary, 
                qtl_group_dir, qtl_report, ref_vcf_dir,
                gwas_col_dict, qtl_col_dict, population, gwas_sample_size, 
                qtl_sample_size, var_id_col_name,
                parallel=False, parallel_worker_num=2, 
                rank_dir=None, current_analysis_order=None, total_numof_analyses=None, 
                whether_schedual=False,
                min_matching_number = 0):
        logging.info(f'Preparing ecaviar gwas files')
        qtl_candidate_df = pd.read_csv(qtl_report, sep=const.column_spliter,
                                        dtype={qtl_col_dict['chrom']: 'category'})
        gwas_candidate_df = pd.read_csv(gwas_cluster_summary, sep=const.column_spliter,
                                        dtype={gwas_col_dict['chrom']: 'category'})
        gwas_cluster_snps_dict = self.__get_cluster_significant_snps_dict(gwas_candidate_df)
        output_base_dir = f'{working_dir}/candidate'
        Path(output_base_dir).mkdir(exist_ok=True, parents=True)

        total_len = len(os.listdir(gwas_cluster_dir))
        ix = 1
        if parallel:
            # Too many workers will cause error "Error: Failed to open"
            with ThreadPoolExecutor(max_workers=parallel_worker_num) as executor:
                futures = []
                for gwas_cluster_file in os.listdir(gwas_cluster_dir):
                    # print(f"gwas_cluster_file: {gwas_cluster_file}")
                    if whether_schedual:
                        outputschedule(rownum=ix,
                            totalnum=total_len,
                            current_analysis_order = current_analysis_order,
                            total_numof_analyses=total_numof_analyses,
                            rank_dir=rank_dir)
                        ix = ix + 1
                    # cluster file name example: chr{chromosome}_{position}-chr{}.tsv.gz
                    result = re.match(r'^(.*)\.tsv(\.gz)?$', gwas_cluster_file)
                    # print(f"ecaviar re result: {result}")
                    if not result or gwas_cluster_file.startswith('.'):
                        continue
                    report_variant_id = result.group(1)
                    # print(f"report_variant_id: {report_variant_id}")
                    # get chr{x}_{position}
                    variant_id = report_variant_id.split('-')[0]
                    # print(f"variant_id: {variant_id}")
                    # gwas文件中显著位点的positions
                    chromosome = utils.get_chromosome_number(report_variant_id.split('-')[1])
                    # print(f"chromosome: {chromosome}")
                    if Path(f'{qtl_group_dir}/{chromosome}').exists():
                        for gene_dir in os.listdir(f'{os.path.join(qtl_group_dir, chromosome)}'):
                            for qtl_file in os.listdir(f"{os.path.join(qtl_group_dir, chromosome, gene_dir)}"):
                                # print(f"qtl_file: {qtl_file}")
                                pheno_id = utils.get_pheno_name(qtl_file)
                                # print(f"pheno_id: {pheno_id}")
                                # qtl一个基因文件中显著位点的positions
                                qtl_significant_positions = self.__get_qtl_significant_snp_positions(
                                    qtl_candidate_df, chromosome, gene_dir, qtl_file)
                                # 如果当前基因文件中至少有一个显著位点存在于当前gwas的显著位点里，才进行gwas和qtl的位点mapping工作
                                if len(set(gwas_cluster_snps_dict[variant_id]) & set(qtl_significant_positions)) < min_matching_number:
                                    continue
                                # if not pheno_id.upper().startswith(const.gene_id_prefix):
                                #     continue
                                gwas_cluster_full_path = f'{gwas_cluster_dir}/{gwas_cluster_file}'
                                qtl_file_full_path = f'{qtl_group_dir}/{chromosome}/{qtl_file}'
                                input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chromosome}.vcf.gz')
                                futures.append(executor.submit(self.process_pheno, working_dir, gwas_cluster_full_path,
                                                            qtl_file_full_path, input_vcf, gwas_col_dict, qtl_col_dict,
                                                            gwas_sample_size, qtl_sample_size, var_id_col_name,
                                                            output_base_dir, pheno_id, variant_id, chromosome, min_matching_number))
                for future in concurrent.futures.as_completed(futures):
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error("".join(traceback.TracebackException.from_exception(exc).format()))
        else:
            for gwas_cluster_file in os.listdir(gwas_cluster_dir):
                # print(f"gwas_cluster_file: {gwas_cluster_file}")
                if whether_schedual:
                    outputschedule(rownum=ix,
                        totalnum=total_len,
                        current_analysis_order = current_analysis_order,
                        total_numof_analyses=total_numof_analyses,
                        rank_dir=rank_dir)
                    ix = ix + 1
                # cluster file name example: chr{chromosome}_{position}-chr{}.tsv.gz
                result = re.match(r'^(.*)\.tsv(\.gz)?$', gwas_cluster_file)
                # print(f"ecaviar re result: {result}")
                if not result or gwas_cluster_file.startswith('.'):
                    continue
                report_variant_id = result.group(1)
                # print(f"report_variant_id: {report_variant_id}")
                # get chr{x}_{position}
                variant_id = report_variant_id.split('-')[0]
                # print(f"variant_id: {variant_id}")
                # gwas文件中显著位点的positions
                chromosome = utils.get_chromosome_number(report_variant_id.split('-')[1])
                # print(f"chromosome: {chromosome}")

                if Path(f'{qtl_group_dir}/{chromosome}').exists():
                    for gene_dir in os.listdir(f'{os.path.join(qtl_group_dir, chromosome)}'):
                        for qtl_file in os.listdir(f"{os.path.join(qtl_group_dir, chromosome, gene_dir)}"):
                            # print(f"qtl_file: {qtl_file}")
                            pheno_id = utils.get_pheno_name(qtl_file)
                            # print(f"pheno_id: {pheno_id}")
                            # qtl一个基因文件中显著位点的positions
                            qtl_significant_positions = self.__get_qtl_significant_snp_positions(
                                qtl_candidate_df, chromosome, gene_dir, qtl_file)
                            # print(f"ecaviar qtl_significant_positions: {qtl_significant_positions}")
                            # 如果当前基因文件中至少有一个显著位点存在于当前gwas的显著位点里，才进行gwas和qtl的位点mapping工作
                            if len(set(gwas_cluster_snps_dict[variant_id]) & set(qtl_significant_positions)) < min_matching_number:
                                print(f"no overlap")
                                continue
                            # if not pheno_id.upper().startswith(const.gene_id_prefix):
                            #     continue
                            gwas_cluster_full_path = f'{gwas_cluster_dir}/{gwas_cluster_file}'
                            # print(f"gwas_cluster_full_path: {gwas_cluster_full_path}")
                            qtl_file_full_path = f'{os.path.join(qtl_group_dir, chromosome, gene_dir, qtl_file)}'
                            # print(f"qtl_file_full_path: {qtl_file_full_path}")
                            input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chromosome}.vcf.gz')
                            self.process_pheno(working_dir, gwas_cluster_full_path, qtl_file_full_path, input_vcf,
                                            gwas_col_dict, qtl_col_dict, gwas_sample_size, qtl_sample_size, var_id_col_name,
                                            output_base_dir, gene_dir, pheno_id, variant_id, chromosome, min_matching_number)
        logging.info(f'Prepare ecaviar input files completed')
        return output_base_dir

    def process_pheno(self, working_dir, gwas_cluster_file_path, qtl_file_path, input_vcf,
                     gwas_col_dict, qtl_col_dict, gwas_sample_size, qtl_sample_size, var_id_col_name,
                     output_base_dir, gene_id, pheno_id, variant_id, chromosome, min_matching_number):
        logging.info(f"sQTL eCAVIAR process phenotype ...")
        gwas_cluster_df = pd.read_csv(gwas_cluster_file_path, sep=const.column_spliter,
                                      usecols=[gwas_col_dict['chrom'], gwas_col_dict['position'],
                                               gwas_col_dict['effect_allele'], gwas_col_dict['other_allele'],
                                               gwas_col_dict['beta'], gwas_col_dict['se'], var_id_col_name],
                                      dtype={
                                          gwas_col_dict['chrom']: 'category',
                                          gwas_col_dict['position']: 'Int64',
                                        #   gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                        #   gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)
                                          gwas_col_dict['effect_allele']: 'category',
                                          gwas_col_dict['other_allele']: 'category'
                                      })
        gwas_cluster_df[self.zscore_col] = \
            gwas_cluster_df[gwas_col_dict['beta']] / gwas_cluster_df[gwas_col_dict['se']]
        # gwas_cluster_df.drop(columns=[gwas_col_dict['beta'], gwas_col_dict['se']], inplace=True)
        qtl_grouped_df = pd.read_csv(qtl_file_path, sep=const.column_spliter,
                                      usecols=[qtl_col_dict['chrom'], 
                                               qtl_col_dict['position'],
                                               qtl_col_dict['alt'], 
                                               qtl_col_dict['ref'],
                                               qtl_col_dict['beta'], 
                                               qtl_col_dict['se'], 
                                               qtl_col_dict['maf'], 
                                               var_id_col_name],
                                      dtype={qtl_col_dict['chrom']: 'category',
                                             qtl_col_dict['position']: 'Int64',
                                            #  qtl_col_dict['alt']: pd.CategoricalDtype(const.SNP_ALLELE),
                                            #  qtl_col_dict['ref']: pd.CategoricalDtype(const.SNP_ALLELE)
                                             qtl_col_dict['alt']: 'category',
                                             qtl_col_dict['ref']: 'category'
                                             })
        # gwas, qtl position matching的交集
        utils.drop_non_intersect_rows(qtl_grouped_df, var_id_col_name, gwas_cluster_df, var_id_col_name)
        if len(gwas_cluster_df) < min_matching_number:
            return
        qtl_grouped_df[self.zscore_col] = \
            qtl_grouped_df[qtl_col_dict['beta']] / qtl_grouped_df[qtl_col_dict['se']]
        # qtl_grouped_df.drop(columns=[qtl_col_dict['beta'], qtl_col_dict['se']], inplace=True)
        vcf_output_dir = os.path.join(working_dir, 'vcf')
        Path(vcf_output_dir).mkdir(parents=True, exist_ok=True)
        candidate_dir = os.path.join(output_base_dir, gene_id, pheno_id, variant_id)
        print(f"candidate_dir: {candidate_dir}")
        Path(candidate_dir).mkdir(parents=True, exist_ok=True)
        # 1.通过gwas, qtl, vcf交集数据生成vcf, run plink 生成LD file
        output_vcf_name = f'{gene_id}_{pheno_id}_{variant_id}.vcf'
        if not Path(input_vcf).exists():
            logging.warning(f'!ref vcf {input_vcf} does not exist')
            return
        gwas_cluster_df.sort_values(by=gwas_col_dict['position'], inplace=True)
        gwas_cluster_df.reset_index(drop=True, inplace=True)
        # utils.adjust_allele_order(gwas_cluster_df,
        #                           gwas_col_dict['effect_allele'],
        #                           gwas_col_dict['other_allele'],
        #                           gwas_col_dict['chrom'],
        #                           gwas_col_dict['position'],
        #                           qtl_grouped_df,
        #                           ref_df_chrom_col_name=qtl_col_dict['chrom'],
        #                           ref_df_pos_col_name=qtl_col_dict['position'],
        #                           ref_df_alt_allele_col_name=qtl_col_dict['alt'],
        #                           ref_df_ref_allele_col_name=qtl_col_dict['ref'],
        #                           gz_col_name=self.zscore_col)
        utils.drop_non_intersect_rows(qtl_grouped_df, var_id_col_name, gwas_cluster_df, var_id_col_name)
        gwas_cluster_df.reset_index(drop=True, inplace=True)
        utils.extract_vcf_data(chromosome, gwas_cluster_df, input_vcf, vcf_output_dir,
                               output_vcf_name, gwas_col_dict['position'], target_snp_col_name=var_id_col_name)
        vcf_matching_file = os.path.join(vcf_output_dir, 'matching', f'{gene_id}_{pheno_id}_{variant_id}.tsv')
        if not os.path.exists(vcf_matching_file) or os.path.getsize(vcf_matching_file) <= 0:
            logging.warning(f'No generated vcf file for gene {pheno_id}')
            return
        vcf_matching_df = pd.read_table(vcf_matching_file, sep=const.column_spliter, header=0,
                                        usecols=[gwas_col_dict['position'], var_id_col_name],
                                        dtype={gwas_col_dict['position']: 'Int64'})

        # vcf_matching_file 包含gwas, qtl, vcf的交集snps
        vcf_matching_df.sort_values(by=gwas_col_dict['position'], inplace=True)
        output_ld_file = f'{candidate_dir}/{gene_id}_{pheno_id}_{variant_id}'
        print(f"sqtl ecaviar output_ld_file: {output_ld_file}")
        os.system(self.shell_command_plink_execute.format(f'{vcf_output_dir}/{output_vcf_name}',
                                                          output_ld_file))
        keeping_col_df = pd.read_table(f'{output_ld_file}.snplist', header=None)
        if keeping_col_df.shape[0] < vcf_matching_df.shape[0]:
            # 从vcf matching snps中移除LD行列为nan的snp
            vcf_matching_df = vcf_matching_df[vcf_matching_df[var_id_col_name].isin(keeping_col_df[0])].copy()
        del keeping_col_df
        # Drop GWAS rows that does not have vcf records
        # SNP in vcf_matching_df is subset of SNP in candidate_gwas_df, so it's fine to drop intersect rows here
        utils.drop_non_intersect_rows(gwas_cluster_df, var_id_col_name, vcf_matching_df, var_id_col_name)
        # Drop qtl rows that does not have vcf records
        utils.drop_non_intersect_rows(qtl_grouped_df, var_id_col_name, vcf_matching_df, var_id_col_name)
        del vcf_matching_df
        # clean vcf file after ld computation
        # utils.delete_file_if_exists(f'{vcf_output_dir}/{output_vcf_name}')

        # 2.根据gwas, qtl, vcf position交集数据生成qtl zscore
        qtl_grouped_df.sort_values(by=qtl_col_dict['position'], inplace=True)
        qtl_grouped_df.reset_index(drop=True, inplace=True)
        qtl_zscore_file = f'{candidate_dir}/qtl_{gene_id}_{pheno_id}_{variant_id}.z'
        # print(f"qtl_zscore_file: {qtl_zscore_file}")
        # qtl_grouped_df.to_csv(qtl_zscore_file, columns=[var_id_col_name, self.zscore_col], header=None,
        #                        sep=' ', index=False)
        qtl_grouped_df['rsid'] = qtl_grouped_df.apply(lambda x: f"{str(x[qtl_col_dict['chrom']])}:{str(x[qtl_col_dict['position']])}", axis=1)
        out_qtl_grouped_df = qtl_grouped_df[['rsid', 
                                         qtl_col_dict['chrom'], 
                                         qtl_col_dict['position'],
                                         qtl_col_dict['alt'], 
                                         qtl_col_dict['ref'],
                                         qtl_col_dict['maf'], 
                                         qtl_col_dict['beta'], 
                                         qtl_col_dict['se']]]
        
        out_qtl_grouped_df.columns = ['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se']
        if len(out_qtl_grouped_df) < min_matching_number:
            logging.info(f"No more than {min_matching_number} SNPs in qtl {gene_id} with lead SNP {variant_id}")
            return
        else:
            out_qtl_grouped_df.to_csv(qtl_zscore_file, header=True, sep=' ', index=False)

        # 3.根据gwas, qtl, vcf position交集数据生成gwas zscore
        del out_qtl_grouped_df
        gwas_zscore_file = f'{candidate_dir}/gwas_{gene_id}_{pheno_id}_{variant_id}.z'
        # print(f"gwas_zscore_file: {gwas_zscore_file}")
        gwas_cluster_df['rsid'] = gwas_cluster_df.apply(lambda x: \
            f"{str(x[gwas_col_dict['chrom']])}:{str(x[gwas_col_dict['position']])}", 
            axis=1)
        gwas_cluster_df = pd.merge(left=gwas_cluster_df,
                        right = qtl_grouped_df[[qtl_col_dict['maf'], var_id_col_name]],
                        left_on = var_id_col_name,
                        right_on = var_id_col_name,
                        how='left')
        out_gwas_cluster_df = gwas_cluster_df[['rsid', 
                                         gwas_col_dict['chrom'], 
                                         gwas_col_dict['position'],
                                         gwas_col_dict['effect_allele'], 
                                         gwas_col_dict['other_allele'],
                                         qtl_col_dict['maf'], 
                                         gwas_col_dict['beta'], 
                                         gwas_col_dict['se']]]
        out_gwas_cluster_df.columns = ['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se']
        if len(out_gwas_cluster_df) < min_matching_number:
            logging.info(f"No more than {min_matching_number} SNPs in gwas {gene_id} with lead SNP {variant_id}")
            return
        else:
            out_gwas_cluster_df.to_csv(gwas_zscore_file, header=True, sep=' ', index=False)
        
        del qtl_grouped_df
        del out_gwas_cluster_df
        del gwas_cluster_df
        # 生成finemap所需的in file
        finemap_in_file = f'{candidate_dir}/{gene_id}_{pheno_id}_{variant_id}.in'
        print(f"finemap_in_file: {finemap_in_file}")
        with open(finemap_in_file, mode='w') as in_file:
            in_file.write('z;ld;snp;config;cred;log;n_samples\n')
            in_file.write(f'{gwas_zscore_file};{output_ld_file}.ld;'
                          f'{candidate_dir}/gwas_{gene_id}_{pheno_id}_{variant_id}.snp;'
                          f'{candidate_dir}/gwas_{gene_id}_{pheno_id}_{variant_id}.config;'
                          f'{candidate_dir}/gwas_{gene_id}_{pheno_id}_{variant_id}.cred;'
                          f'{candidate_dir}/gwas_log.log;{gwas_sample_size}\n')
            in_file.write(f'{qtl_zscore_file};{output_ld_file}.ld;'
                          f'{candidate_dir}/qtl_{gene_id}_{pheno_id}_{variant_id}.snp;'
                          f'{candidate_dir}/qtl_{gene_id}_{pheno_id}_{variant_id}.config;'
                          f'{candidate_dir}/qtl_{gene_id}_{pheno_id}_{variant_id}.cred;'
                          f'{candidate_dir}/qtl_log.log;{qtl_sample_size}')

    def __get_cluster_significant_snps_dict(self, cluster_df):
        cluster_snps_dict = {}
        for _, row in cluster_df.iterrows():
            cluster_snps_dict[row.loc['lead_SNP']] = self.__convert_positions_str_to_list(row.loc['positions'])
        return cluster_snps_dict

    def __get_qtl_significant_snp_positions(self, qtl_df, chrom, gene_id, gene_file_name):
        results = qtl_df[(qtl_df['pheno_file'] == gene_file_name) & \
                          (qtl_df['gene'] == gene_id) & \
                          (qtl_df['chrom'].astype(str) == str(chrom))]
        if results.empty:
            return []
        else:
            return self.__convert_positions_str_to_list(results.iloc[0]['positions'])

    def __convert_positions_str_to_list(self, positions_str):
        if isinstance(positions_str, str):
            return json.loads(positions_str)
        elif isinstance(positions_str, list):
            return positions_str
        else:
            return []


if __name__ == '__main__':
    cfg = config.ConfigHolder()
    glob_processor = gdp.Processor(cfg_holder=cfg)
    _working_dir = os.path.join('/Users/nicklin/Desktop/CAD/test', ECaviarDataProcessor.ECAVIAR_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    pop = 'EUR'
    gwas_data_processor = ECaviarDataProcessor()

    _gwas_col_dict = {'snp': 'oldID', 'chrom': 'chr', 'position': 'hm_pos', 'beta': 'beta',
                      'effect_allele': 'a1', 'other_allele': 'a2', 'pvalue': 'pval',
                      'se': 'se'}
    _qtl_col_dict = {'snp': 'rsid', 'chrom': 'chromosome', 'position': 'position', 'beta': 'beta', 'alt': 'alt',
                      'ref': 'ref', 'pvalue': 'pvalue', 'se': 'se', 'pheno_id': 'molecular_trait_id', 'maf': 'maf'}

    gwas_data_processor.prepare(working_dir=_working_dir,
                                gwas_cluster_dir='/Volumes/HD/biodata/colocalization-tools/preprocessed/gwas/cad/clustered',
                                gwas_cluster_summary='/Volumes/HD/biodata/colocalization-tools/preprocessed/gwas/cad/cluster_summary.tsv',
                                qtl_group_dir='/Volumes/HD/biodata/colocalization-tools/preprocessed/qtl/Artery_Aorta/grouped',
                                qtl_report='/Volumes/HD/biodata/colocalization-tools/preprocessed/qtl/Artery_Aorta/filtered_gene.tsv',
                                ref_vcf_dir='/Volumes/HD/biodata/colocalization-tools/raw/vcf/hg38_phased',
                                gwas_col_dict=_gwas_col_dict,
                                qtl_col_dict=_qtl_col_dict,
                                population=pop,
                                gwas_sample_size=296525,
                                qtl_sample_size=948,
                                var_id_col_name='var_id')
