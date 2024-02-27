import concurrent
import logging
import os
import re
import shutil
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import coloc_utils as utils, global_data_process as gdp, constants as const


class TWAS:
    COLOC_TOOL_NAME = 'twas'

    def __init__(self):
        logging.info('init TWAS')

    def run(self,
            working_dir=None,
            weights_path=None,
            gwas_chrom_group_dir=None,
            gwas_col_dict=None,
            ref_vcf_dir=None,
            population=None,
            parallel=False,
            tools_config_file=None,
            parallel_worker_num=4):
        
        
        start_time = datetime.now()
        logging.info(f'run_twas start at: {start_time}')
        output_file = self.__get_output_file(working_dir)
        Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
        input_dir = os.path.join(working_dir, 'input')
        Path(input_dir).mkdir(parents=True, exist_ok=True)
        twas_path = shutil.which("FUSION.assoc_test.R")
        if twas_path is None:
            raise ValueError(f'FUSION.assoc_test.R not found, '
                             f'did you config it in PATH env(or did you activate conda env?)')
        if weights_path is None:
            raise ValueError(f'TWAS weight files(.pos file and .RDat files) must be provided to run TWAS')
        ld_ref_prefix = os.path.join(ref_vcf_dir, population.upper(), f'chr')
        if parallel:
            with ThreadPoolExecutor(max_workers=parallel_worker_num) as executor:
                futures = []
                for chrom_group_file in os.listdir(gwas_chrom_group_dir):
                    chrom = self.__ensure_chrom_ld_ref(chrom_group_file, ref_vcf_dir, population)
                    if chrom is None:
                        continue
                    chrom_group_file_path = os.path.join(gwas_chrom_group_dir, chrom_group_file)
                    futures.append(
                        executor.submit(self.process_chrom, chrom, chrom_group_file_path, ld_ref_prefix,
                                        gwas_col_dict, input_dir, twas_path, weights_path, tools_config_file))
                for future in concurrent.futures.as_completed(futures):
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error("".join(traceback.TracebackException.from_exception(exc).format()))

        else:
            for chrom_group_file in os.listdir(gwas_chrom_group_dir):
                chrom = self.__ensure_chrom_ld_ref(chrom_group_file, ref_vcf_dir, population)
                if chrom is None:
                    continue
                chrom_group_file_path = os.path.join(gwas_chrom_group_dir, chrom_group_file)
                self.process_chrom(chrom, chrom_group_file_path, ld_ref_prefix, gwas_col_dict, input_dir,
                                   twas_path, weights_path, tools_config_file)

        self.__analyze_result(input_dir, output_file)
        if not os.path.exists(output_file) or os.path.getsize(output_file) <= 0:
            logging.warning(f'Process completed, duration {datetime.now() - start_time}, no result found')
        else:
            logging.info(
                f'Process completed, duration {datetime.now() - start_time}, check {output_file} for result!')
        return output_file

    def __ensure_chrom_ld_ref(self, chrom_group_file, ref_vcf_dir, population):
        
        
        if not (chrom_group_file.endswith('.tsv') or chrom_group_file.endswith('.tsv.gz')):
            return None
        chr_nums = re.findall(r'\d+', chrom_group_file)
        if len(chr_nums) == 0:
            return None
        chrom = chr_nums[0]
        chrom_ld_ref_prefix = os.path.join(ref_vcf_dir, population.upper(), f'chr{chrom}')
        if not os.path.exists(f'{chrom_ld_ref_prefix}.bed') or not os.path.exists(
                f'{chrom_ld_ref_prefix}.bim') or not os.path.exists(f'{chrom_ld_ref_prefix}.fam'):
            logging.info(f'plink binary ref LD files for chrom {chrom} do not exist, trying to generate from vcf')
            input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chrom}.vcf.gz')
            if not os.path.exists(input_vcf):
                logging.warning(f'!ref vcf {input_vcf} does not exist')
                return None
            # https://github.com/gusevlab/fusion_twas/issues/28
            os.system(f'plink --silent --vcf {input_vcf} '
                      f'--maf 0.01 --geno 0.01 --mind 0.01 --hwe 1e-06 --snps-only just-acgt '
                      f'--biallelic-only strict --make-bed --out {chrom_ld_ref_prefix}')
        return chrom

    def process_chrom(self, chrom, gwas_file, ld_ref_prefix, gwas_col_dict, input_dir, twas_path,
                      weights_path, tools_config_file):
        
        
        logging.warning(f'Processing for chromosome {chrom} start')
        gwas_chrom_df = pd.read_table(gwas_file, sep=const.column_spliter,
                                      usecols=[gwas_col_dict['snp'], gwas_col_dict['chrom'], gwas_col_dict['position'],
                                               gwas_col_dict['effect_allele'], gwas_col_dict['other_allele'],
                                               gwas_col_dict['beta'], gwas_col_dict['se']],
                                      dtype={gwas_col_dict['chrom']: 'category',
                                             gwas_col_dict['position']: 'Int64',
                                             gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                             gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE)})
        if gwas_chrom_df.empty:
            logging.warning(f'gwas input size is 0')
            return
        gwas_chrom_df['Z'] = gwas_chrom_df[gwas_col_dict['beta']] / gwas_chrom_df[gwas_col_dict['se']]
        gwas_chrom_df.drop(columns=[gwas_col_dict['chrom'], gwas_col_dict['position'],
                                    gwas_col_dict['beta'], gwas_col_dict['se']], inplace=True)
        #  NOTE!! TWAS input must contain SNP column,
        #  and values in this column must match values in the 2nd column of .bim file (indicated by --ref_ld_chr)
        gwas_chrom_df.rename({gwas_col_dict['snp']: 'SNP', gwas_col_dict['effect_allele']: 'A1',
                              gwas_col_dict['other_allele']: 'A2'}, axis='columns', inplace=True)
        gwas_chrom_df = gwas_chrom_df.reindex(columns=['SNP', 'A1', 'A2', 'Z'], copy=False)
        gwas_chrom_input = os.path.join(input_dir, f'chr{chrom}.tsv.gz')
        gwas_chrom_df.to_csv(gwas_chrom_input, sep=const.column_spliter, header=True, index=False)
        del gwas_chrom_df
        chrom_out_file = os.path.join(input_dir, f'chr{chrom}_result.tsv')
        custom_params = utils.get_tools_params(self.COLOC_TOOL_NAME, tools_config_file=tools_config_file)
        if custom_params is None:
            custom_params = ''
        # Error "cannot open file '/xxx/yyy/chr.bim': No such file or directory" can be ignored" can be ignored
        # it means no weight data of this chrom is available in .pos file
        os.system(f'Rscript --no-save --no-restore {twas_path} '
                  f'--sumstats  {gwas_chrom_input} '
                  f'--weights {weights_path} '
                  f'--weights_dir {os.path.dirname(weights_path)} '
                  f'--ref_ld_chr {ld_ref_prefix} '
                  f'--chr {chrom} '
                  f'--max_impute 1 '
                  f'--out {chrom_out_file} '
                  f'{custom_params}')
        logging.warning(f'Processing for chromosome {chrom} completed')

    def __get_output_file(self, working_dir):
        
        
        _output_file_name = f'{TWAS.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        return os.path.join(working_dir, 'analyzed', _output_file_name)

    def __analyze_result(self, output_dir,
                         final_report_file):
        # Merge every single result
        
        
        single_result_list = []
        for single_result in os.listdir(output_dir):
            if not single_result.endswith('_result.tsv'):
                continue
            single_result_list.append(pd.read_table(os.path.join(output_dir, single_result)))
        if len(single_result_list) == 0:
            return
        report_df = pd.concat(single_result_list)
        if report_df.shape[0] == 0:
            return
        report_df.replace(to_replace=r'^\s*NA\s*$', value=pd.NA, inplace=True, regex=True)
        report_df.dropna(subset='TWAS.P', inplace=True)
        if report_df.shape[0] == 0:
            return
        report_df['TWAS.P'] = report_df['TWAS.P'].astype(float)
        gene_df = report_df['ID'].str.split('.', n=2, expand=True)
        report_df['ID'] = gene_df[0]
        report_df.sort_values(by='TWAS.P', ascending=True, inplace=True)
        report_df.rename(columns={'ID': 'gene_id'}, inplace=True)
        report_df.drop_duplicates(subset=['gene_id'], inplace=True)
        if len(report_df) <= 0:
            utils.delete_file_if_exists(final_report_file)
        else:
            report_df['TWAS.P'] = report_df['TWAS.P'].apply(lambda x: "{:.3e}".format(x))
            report_df.to_csv(final_report_file, sep=const.output_spliter, header=True, index=False)


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    _gwas_chrom_group_dir = glob_processor.gwas_output_dir
    if len(os.listdir(_gwas_chrom_group_dir)) == 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    weight_pos_file = utils.get_twas_ref_files(glob_processor.global_config)
    if not os.path.exists(weight_pos_file) or os.path.getsize(weight_pos_file) <= 0:
        raise ValueError(f'TWAS weight pos file not found')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, TWAS.COLOC_TOOL_NAME)
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    twas = TWAS()
    twas.run(_working_dir,
             weight_pos_file,
             _gwas_chrom_group_dir,
             glob_processor.gwas_col_dict,
             glob_processor.ref_vcf_dir,
             pop,
             tools_config_file='/Users/nicklin/gitrepo/bio/colocalization-tools/resource/tools_config.yml',
             parallel=True)
