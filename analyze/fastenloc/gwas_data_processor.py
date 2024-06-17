import logging
import os
from datetime import datetime
from pathlib import Path
import sys
import pandas as pd

from common import global_data_process as gdp, constants as const, coloc_utils as utils

def outputschedule(rownum, numofeqtlloci,currenttissuenum, numoftissues, rank_dir):
    calculated_schedule = int(rownum/numofeqtlloci * 80/numoftissues + 80/numoftissues * (currenttissuenum - 1))
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()


class FastenlocGwasProcessor:

    def __init__(self):
        logging.info('init FastenlocGwasProcessor')

    def __loc_array_init(self, chr_str, chr_len):
        arr = [f'loc_{chr_str}_{x}' for x in range(chr_len)]
        return arr

    def get_output_gwas_dir(self, output_base_dir):
        return f'{output_base_dir}/gwas'

    def get_output_torus_output(self, output_gwas_dir):
        return f'{output_gwas_dir}/pip'

    def get_output_torus_output_file(self, output_torus_output_dir):
        return f'{output_torus_output_dir}/torus_output.pip'

    def prepare_gwas_data(self, working_dir=None, gwas_preprocessed_file=None,
                          gwas_col_dict=None, ld_block_loci_file=None, rank_dir=None,
                          currenttissuenum=None, numoftissues=None, 
                          whether_schedual=False):
        start_time = datetime.now()
        logging.info(f'Grouping fastenloc preprocessed gwas file at {start_time}')

        output_base_dir = working_dir
        output_gwas_dir = self.get_output_gwas_dir(output_base_dir)
        output_torus_input_dir = f'{output_gwas_dir}/torus_input'
        output_torus_input_file = f'{output_torus_input_dir}/torus_input'
        utils.delete_file_if_exists(output_torus_input_file)
        output_torus_output_dir = self.get_output_torus_output(output_gwas_dir)
        output_torus_output_file = self.get_output_torus_output_file(output_torus_output_dir)

        dap_use_col = [gwas_col_dict['beta'],
                       gwas_col_dict['se'],
                       gwas_col_dict['position'], gwas_col_dict['chrom'],
                       gwas_col_dict['effect_allele'],
                       gwas_col_dict['other_allele'], 'ref_', 'alt_']

        shell_command_torus_execute = 'torus -d {} --load_zval -dump_pip {}'
        shell_compress_file = 'gzip -k -f {}'

        Path(output_torus_input_dir).mkdir(parents=True, exist_ok=True)
        Path(output_torus_output_dir).mkdir(parents=True, exist_ok=True)
        gwas_snp_count = 0
        tmp = pd.read_csv(gwas_preprocessed_file, sep=const.column_spliter)
        total_num = len(tmp)
        del tmp
        ix = 0
        with pd.read_table(gwas_preprocessed_file, sep=const.column_spliter, usecols=dap_use_col,
                           dtype={gwas_col_dict['position']: 'Int64',
                                  gwas_col_dict['chrom']: 'category',
                                  gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                  gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
                                  'ref_': pd.CategoricalDtype(const.SNP_ALLELE),
                                  'alt_': pd.CategoricalDtype(const.SNP_ALLELE)},
                           iterator=True, chunksize=500000) as reader:
            for chunk in reader:
                ix = ix + len(chunk)
                if whether_schedual == True:
                    outputschedule(rownum=ix,
                                numofeqtlloci=total_num,
                                currenttissuenum = currenttissuenum,
                                numoftissues=numoftissues,
                                rank_dir=rank_dir)
                gwas_snp_count += chunk.shape[0]
                chunk['zscore'] = chunk[gwas_col_dict['beta']] / chunk[gwas_col_dict['se']]
                chunk.drop(columns=[gwas_col_dict['beta'], gwas_col_dict['se']], inplace=True)
                # Adjust allele order
                ref_df = chunk[[gwas_col_dict['chrom'], gwas_col_dict['position'], 'alt_', 'ref_']]
                logging.info(f'Aligning allele order')
                utils.adjust_allele_order(chunk,
                                          gwas_col_dict['effect_allele'],
                                          gwas_col_dict['other_allele'],
                                          gwas_col_dict['chrom'],
                                          gwas_col_dict['position'],
                                          ref_df,
                                          ref_df_chrom_col_name=gwas_col_dict['chrom'],
                                          ref_df_pos_col_name=gwas_col_dict['position'],
                                          ref_df_alt_allele_col_name='alt_',
                                          ref_df_ref_allele_col_name='ref_',
                                          gz_col_name='zscore',
                                          drop_ref_df_non_intersect_items=False)
                logging.info(f'Align allele order completed')
                del ref_df
                # chr1_13550_G_A_b38
                chunk['variant_id'] = 'chr' + chunk[gwas_col_dict['chrom']].astype(str) + '_' + \
                                      chunk[gwas_col_dict['position']].astype(str) + '_' + \
                                      chunk['ref_'].astype(str) + '_' + chunk['alt_'].astype(str)
                chunk.drop(columns=[gwas_col_dict['effect_allele'], gwas_col_dict['other_allele'], 'ref_', 'alt_'],
                           inplace=True)
                logging.info(f'Adding zscore/variant_id completed, memory usage: {chunk.memory_usage(deep=True)}')
                # get loc range in eur_ld.hg38.bed for loc column
                df_ld_bed_loc = pd.read_csv(ld_block_loci_file, sep=const.column_spliter,
                                            dtype={'start': 'Int64', 'stop': 'Int64'})
                df_ld_bed_loc.chr = df_ld_bed_loc.chr.apply(lambda x: int(x.replace('chr', '')))

                df_gwas_chr_list = list(chunk.groupby(gwas_col_dict['chrom']))
                df_gwas_list = []
                for gwas_chr in df_gwas_chr_list:
                    chr_num = int(gwas_chr[0])
                    gwas_chr_df = gwas_chr[1]
                    bins_list = df_ld_bed_loc[df_ld_bed_loc['chr'] == chr_num].loc[:, 'start'] - 1
                    gwas_chr_df['loc'] = pd.cut(x=gwas_chr_df[gwas_col_dict['position']],
                                                bins=list(bins_list),
                                                labels=self.__loc_array_init(chr_num, len(bins_list) - 1))
                    gwas_chr_df.drop(columns=[gwas_col_dict['chrom'], gwas_col_dict['position']], inplace=True)
                    df_gwas_list.append(gwas_chr_df)
                if os.path.exists(output_torus_input_file) and os.path.getsize(output_torus_input_file) > 0:
                    mode = 'a'
                    header = False
                else:
                    mode = 'w'
                    header = True
                df_input_torus = pd.concat(df_gwas_list)
                # torus input file include columns variant_id、loc、zscore
                df_input_torus.to_csv(output_torus_input_file, columns=['variant_id', 'loc', 'zscore'],
                                      sep=const.output_spliter, mode=mode, header=header, index=False)
                del df_gwas_list, df_input_torus
        logging.info(f'prepare fastenloc gwas file all columns duration {datetime.now() - start_time}')
        os.system(shell_compress_file.format(output_torus_input_file))
        #
        os.system(shell_command_torus_execute.format(f'{output_torus_input_file}.gz', output_torus_output_file))
        os.system(shell_compress_file.format(output_torus_output_file))
        logging.info(f'prepare fastenloc gwas file at: {datetime.now()}, duration： {datetime.now() - start_time}')
        return output_torus_output_file, gwas_snp_count
        # clean temp file folder
        # utils.delete_dir(output_torus_input_dir)


if __name__ == '__main__':
    fastenloc_gwas_pro = FastenlocGwasProcessor()
    processor = gdp.Processor()
    _working_dir = os.path.join(processor.tool_parent_dir, 'fastenloc')
    fastenloc_gwas_pro.prepare_gwas_data(working_dir=_working_dir,
                                         gwas_preprocessed_file=processor.gwas_preprocessed_file,
                                         gwas_col_dict=processor.gwas_col_dict,
                                         ld_block_loci_file=processor.global_config['input']['ld_block_loci_file'])
