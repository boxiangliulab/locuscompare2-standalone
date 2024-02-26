import os
import sys
from pathlib import Path
from datetime import datetime
import logging
import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import coloc_utils as utils, global_data_process as gdp, constants as const


class PredixcanGwasProcessor:
    COLOC_TOOL_NAME = 'predixcan'
    PREDIXCAN_VAR_ID_COL_NAME = 'p_var_id'

    def __init__(self):
        logging.info('init PredixcanGwasProcessor')

    def prepare(self,
                working_dir=None,
                gwas_preprocessed_file=None,
                var_id_col_name=None,
                gwas_col_dict=None):
        
        
        start_time = datetime.now()
        logging.info(f'Preparing gwas file at {start_time}')
        output_processed_file = self.__get_output_file(working_dir)
        utils.delete_file_if_exists(output_processed_file)
        with pd.read_table(gwas_preprocessed_file, sep=const.column_spliter,
                           usecols=[gwas_col_dict['chrom'], gwas_col_dict['position'], gwas_col_dict['effect_allele'],
                                    gwas_col_dict['other_allele'], gwas_col_dict['beta'], gwas_col_dict['se'], 'alt_',
                                    'ref_'],
                           dtype={gwas_col_dict['position']: 'Int64',
                                  gwas_col_dict['chrom']: 'category',
                                  gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                  gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                  'ref_': pd.CategoricalDtype(const.SNP_ALLELE),
                                  'alt_': pd.CategoricalDtype(const.SNP_ALLELE)},
                           iterator=True, chunksize=500000, header=0) as reader:
            for chunk in reader:
                chunk[PredixcanGwasProcessor.PREDIXCAN_VAR_ID_COL_NAME] = \
                    'chr' + chunk[gwas_col_dict['chrom']].astype(str) \
                    + '_' + chunk[gwas_col_dict['position']].astype(str) \
                    + '_' + chunk['ref_'].astype(str) \
                    + '_' + chunk['alt_'].astype(str) \
                    + '_b38'
                chunk.drop(columns=[gwas_col_dict['chrom'], gwas_col_dict['position'], 'alt_', 'ref_'], inplace=True)
                if os.path.exists(output_processed_file) and os.path.getsize(output_processed_file) > 0:
                    mode = 'a'
                    header = False
                else:
                    mode = 'w'
                    header = True
                chunk[[PredixcanGwasProcessor.PREDIXCAN_VAR_ID_COL_NAME,
                       gwas_col_dict['effect_allele'],
                       gwas_col_dict['other_allele'],
                       gwas_col_dict['beta'],
                       gwas_col_dict['se']]
                ].to_csv(output_processed_file, mode=mode, sep=const.output_spliter, header=header, index=False)
        logging.info(f'Preparing gwas file {output_processed_file} completed at {datetime.now()},'
                     f'duration {datetime.now() - start_time}')
        return output_processed_file

    def __get_output_file(self, working_dir):
        
        
        return os.path.join(working_dir, 'gwas_processed.tsv.gz')


if __name__ == '__main__':
    glob_processor = gdp.Processor()
    if not os.path.exists(glob_processor.gwas_preprocessed_file) or os.path.getsize(
            glob_processor.gwas_preprocessed_file) <= 0:
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    _working_dir = os.path.join(glob_processor.tool_parent_dir, PredixcanGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    gwas_data_processor = PredixcanGwasProcessor()
    gwas_data_processor.prepare(_working_dir,
                                glob_processor.gwas_preprocessed_file,
                                gdp.Processor.VAR_ID_COL_NAME,
                                glob_processor.gwas_col_dict)
