import os
import sys
from pathlib import Path
from datetime import datetime
import logging
import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import coloc_utils as utils, global_data_process as gdp, constants as const


class SmrGwasProcessor:
    COLOC_TOOL_NAME = 'smr'

    def __init__(self):
        logging.info('init SmrGwasProcessor')

    def prepare(self,
                working_dir=None,
                gwas_preprocessed_file=None,
                gwas_col_dict=None):
        start_time = datetime.now()
        logging.info(f'Grouping preprocessed gwas file by chrom at {start_time}')
        # Remove files if they already exist
        gwas_chrom_group_dir = self.__get_output_dir(working_dir)
        Path(gwas_chrom_group_dir).mkdir(exist_ok=True, parents=True)
        gwas_chrom_group_file_pattern = self.__get_output_file_pattern()
        for file in os.listdir(gwas_chrom_group_dir):
            utils.delete_file_if_exists(os.path.join(gwas_chrom_group_dir, file))
        # Process by batch
        with pd.read_table(gwas_preprocessed_file, sep=const.column_spliter,
                           dtype={gwas_col_dict['chrom']: 'category',
                                  gwas_col_dict['position']: 'Int64'},
                           iterator=True, chunksize=500000, header=0) as reader:
            for chunk in reader:
                grouped = chunk.groupby(gwas_col_dict['chrom'])
                for name, group in grouped:
                    group_file = os.path.join(gwas_chrom_group_dir,
                                              gwas_chrom_group_file_pattern.format(name))
                    if os.path.exists(group_file) and os.path.getsize(group_file) > 0:
                        mode = 'a'
                        header = False
                    else:
                        mode = 'w'
                        header = True
                    group.to_csv(group_file, mode=mode, sep=const.output_spliter, header=header, index=False)
        logging.info(f'Grouping preprocessed gwas file to {gwas_chrom_group_dir} by chrom completed at {datetime.now()},'
              f'duration {datetime.now() - start_time}')
        return self.__get_output_dir(working_dir), self.__get_output_file_pattern()

    def __get_output_dir(self, working_dir):
        return os.path.join(working_dir, 'gwas_grouped')

    def __get_output_file_pattern(self):
        return 'chr{}.tsv'


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    if not os.path.exists(glob_processor.gwas_preprocessed_file) or os.path.getsize(
            glob_processor.gwas_preprocessed_file) <= 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, SmrGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    smr_gwas_processor = SmrGwasProcessor()
    smr_gwas_processor.prepare(_working_dir,
                               glob_processor.gwas_preprocessed_file,
                               glob_processor.gwas_col_dict)
