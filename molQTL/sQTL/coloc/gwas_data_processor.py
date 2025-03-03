import os
import sys
from pathlib import Path
from datetime import datetime
import logging
import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import utils, global_data_process as gdp, constants as const


class ColocGwasProcessor:
    COLOC_TOOL_NAME = 'coloc'

    def __init__(self):
        logging.info('init sQTL ColocGwasProcessor')

    def prepare(self,
                working_dir=None,
                gwas_cluster_output_dir=None,
                population=None,
                ref_vcf_dir=None,
                var_id_col_name=None,
                gwas_col_dict=None):
        
        
        start_time = datetime.now()
        vcf_output_dir = self.get_output_vcf_dir(working_dir)
        utils.delete_dir(vcf_output_dir)
        Path(vcf_output_dir).mkdir(exist_ok=True, parents=True)
        logging.info(f'start to extract vcf data at {start_time}')

        for cluster_file_name in os.listdir(gwas_cluster_output_dir):
            cluster_file_path = f'{gwas_cluster_output_dir}/{cluster_file_name}'
            print(f"cluster_file_path: {cluster_file_path}")
            ref_snp_rsid = cluster_file_name.split('-')[0]
            chromosome = cluster_file_name.split('-')[1].split('.')[0]
            output_name = self.get_output_vcf_file_pattern().format(ref_snp_rsid, chromosome)
            input_vcf = os.path.join(ref_vcf_dir, population.upper(), f'chr{chromosome}.vcf.gz')
            if Path(input_vcf).exists():
                snp_positions_rsids = pd.read_csv(cluster_file_path, sep=const.column_spliter,
                                                  usecols=[var_id_col_name, gwas_col_dict['position']])
                utils.extract_vcf_data(chromosome, snp_positions_rsids, input_vcf, vcf_output_dir,
                                       output_name,
                                       gwas_col_dict['position'],
                                       target_snp_col_name=var_id_col_name)
            else:
                logging.warning(f'!ref vcf {input_vcf} does not exist')
        logging.info(f'Extract vcf for gwas clustered file completed at {datetime.now()}, '
              f'destination directory {vcf_output_dir}'
              f'duration {datetime.now() - start_time}')
        return vcf_output_dir

    def get_output_vcf_dir(self, working_dir):
        
        
        return os.path.join(working_dir, 'vcf')

    def get_output_vcf_file_pattern(self):
        
        
        return '{}-{}.vcf'


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    if len(os.listdir(glob_processor.gwas_cluster_output_dir)) == 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, ColocGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    gwas_data_processor = ColocGwasProcessor()
    gwas_data_processor.prepare(_working_dir,
                                glob_processor.gwas_cluster_output_dir,
                                pop,
                                glob_processor.ref_vcf_dir,
                                gdp.Processor.VAR_ID_COL_NAME,
                                glob_processor.gwas_col_dict)
