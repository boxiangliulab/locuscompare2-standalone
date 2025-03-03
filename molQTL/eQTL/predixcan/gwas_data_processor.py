import os
import sys
from pathlib import Path
from datetime import datetime
import logging
import pandas as pd

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import utils, global_data_process as gdp, constants as const

def outputschedule(rownum, totalnum, current_analysis_order, total_numof_analyses, rank_dir):
    calculated_schedule = int(rownum/totalnum * 80/total_numof_analyses + 80/total_numof_analyses * (current_analysis_order - 1))
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()

class PredixcanGwasProcessor:
    COLOC_TOOL_NAME = 'predixcan'
    PREDIXCAN_VAR_ID_COL_NAME = 'p_var_id'

    def __init__(self):
        logging.info('init PredixcanGwasProcessor')

    def prepare(self,
                working_dir=None,
                gwas_preprocessed_file=None,
                var_id_col_name=None,
                gwas_col_dict=None, rank_dir=None,
                current_analysis_order=None, total_numof_analyses=None, 
                whether_schedual=False):
        
        start_time = datetime.now()
        logging.info(f'Preparing gwas file at {start_time}')
        output_processed_file = self.__get_output_file(working_dir)
        utils.delete_file_if_exists(output_processed_file)
        if Path(gwas_preprocessed_file).is_dir():
            ## each gwas loci
            loci = os.listdir(gwas_preprocessed_file)
            total_num = len(loci)
            ix = 0
            output_processed_file_ls = []
            for eachgwasloci in loci:
                output_processed_file = self.__get_loci_output_file(working_dir, eachgwasloci)
                utils.delete_file_if_exists(output_processed_file)
                with pd.read_table(os.path.join(gwas_preprocessed_file, eachgwasloci), sep=const.column_spliter,
                                usecols=[gwas_col_dict['chrom'], gwas_col_dict['position'], gwas_col_dict['effect_allele'],
                                            gwas_col_dict['other_allele'], gwas_col_dict['beta'], gwas_col_dict['se'], 'alt_',
                                            'ref_'],
                                dtype={gwas_col_dict['position']: 'Int64',
                                        gwas_col_dict['chrom']: 'category',
                                        #   gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                        #   gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                        #   'ref_': pd.CategoricalDtype(const.SNP_ALLELE),
                                        #   'alt_': pd.CategoricalDtype(const.SNP_ALLELE)
                                        gwas_col_dict['effect_allele']: 'category',
                                        gwas_col_dict['other_allele']: 'category',
                                        'ref_': 'category',
                                        'alt_': 'category'
                                        },
                                iterator=True, chunksize=500000, header=0) as reader:
                    for chunk in reader:
                        ix = ix + len(chunk)
                        if whether_schedual == True:
                            outputschedule(rownum=ix,
                                        totalnum=total_num,
                                        current_analysis_order = current_analysis_order,
                                        total_numof_analyses=total_numof_analyses,
                                        rank_dir=rank_dir)
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
                output_processed_file_ls.append(output_processed_file)
            return output_processed_file_ls
        else:
            tmp = pd.read_csv(gwas_preprocessed_file, sep=const.column_spliter)
            total_num = len(tmp)
            del tmp
            ix = 0
            with pd.read_table(gwas_preprocessed_file, sep=const.column_spliter,
                            usecols=[gwas_col_dict['chrom'], gwas_col_dict['position'], gwas_col_dict['effect_allele'],
                                        gwas_col_dict['other_allele'], gwas_col_dict['beta'], gwas_col_dict['se'], 'alt_',
                                        'ref_'],
                            dtype={gwas_col_dict['position']: 'Int64',
                                    gwas_col_dict['chrom']: 'category',
                                    #   gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                    #   gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                    #   'ref_': pd.CategoricalDtype(const.SNP_ALLELE),
                                    #   'alt_': pd.CategoricalDtype(const.SNP_ALLELE)
                                    gwas_col_dict['effect_allele']: 'category',
                                    gwas_col_dict['other_allele']: 'category',
                                    'ref_': 'category',
                                    'alt_': 'category'
                                    },
                            iterator=True, chunksize=500000, header=0) as reader:
                for chunk in reader:
                    ix = ix + len(chunk)
                    if whether_schedual == True:
                        outputschedule(rownum=ix,
                                    totalnum=total_num,
                                    current_analysis_order = current_analysis_order,
                                    total_numof_analyses=total_numof_analyses,
                                    rank_dir=rank_dir)
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

    def __get_loci_output_file(self, working_dir, eachgwasloci):
        return os.path.join(working_dir, eachgwasloci)
    
    
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
