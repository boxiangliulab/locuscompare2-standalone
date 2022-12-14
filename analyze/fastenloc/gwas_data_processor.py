import pandas as pd
from common import global_data_process as gdp, constants as const
import os
from pathlib import Path
from datetime import datetime
import logging


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
        return f'{output_torus_output_dir}/tours_output.pip'

    def prepare_gwas_data(self, working_dir=None, gwas_preprocessed_file=None,
                          gwas_col_dict=None, var_id_col_name=None, ld_block_loci_file=None):
        start_time = datetime.now()
        logging.info(f'Grouping fastenloc preprocessed gwas file at {start_time}')

        output_base_dir = working_dir
        output_gwas_dir = self.get_output_gwas_dir(output_base_dir)
        output_torus_input_dir = f'{output_gwas_dir}/torus_input'
        output_torus_input_file = f'{output_torus_input_dir}/tours_input'
        output_torus_output_dir = self.get_output_torus_output(output_gwas_dir)
        output_torus_output_file = self.get_output_torus_output_file(output_torus_output_dir)

        dap_use_col = [var_id_col_name, gwas_col_dict['beta'],
                       gwas_col_dict['se'],
                       gwas_col_dict['position'], gwas_col_dict['chrom'],
                       gwas_col_dict['effect_allele'],
                       gwas_col_dict['other_allele']]

        shell_command_torus_execute = 'torus -d {} --load_zval -dump_pip {}'
        shell_compress_file = 'gzip -k -f {}'

        Path(output_torus_input_dir).mkdir(parents=True, exist_ok=True)
        Path(output_torus_output_dir).mkdir(parents=True, exist_ok=True)

        df_gwas_file = pd.read_table(gwas_preprocessed_file, sep=const.column_spliter, usecols=dap_use_col)

        df_gwas_file['zscore'] = df_gwas_file[gwas_col_dict['beta']] / df_gwas_file[gwas_col_dict['se']]
        # chr1_13550_G_A_b38
        df_gwas_file['variant_id'] = 'chr' + df_gwas_file[gwas_col_dict['chrom']].astype(str) + '_' + \
                                     df_gwas_file[gwas_col_dict['position']].map(
                                         str) + '_' + df_gwas_file[gwas_col_dict['other_allele']] + '_' + \
                                     df_gwas_file[gwas_col_dict['effect_allele']]

        # get loc range in eur_ld.hg38.bed for loc column
        df_ld_bed_loc = pd.read_csv(ld_block_loci_file, sep=const.column_spliter,
                                    dtype={'start': 'Int64', 'stop': 'Int64'})
        df_ld_bed_loc.chr = df_ld_bed_loc.chr.apply(lambda x: int(x.replace('chr', '')))

        df_gwas_chr_list = list(df_gwas_file.groupby(gwas_col_dict['chrom']))
        df_gwas_list = []
        for gwas_chr in df_gwas_chr_list:
            chr_num = int(gwas_chr[0])
            gwas_chr_df = gwas_chr[1]
            bins_list = df_ld_bed_loc[df_ld_bed_loc['chr'] == chr_num].loc[:, 'start'] - 1
            gwas_chr_df['loc'] = pd.cut(x=gwas_chr_df[gwas_col_dict['position']],
                                        bins=list(bins_list),
                                        labels=self.__loc_array_init(chr_num, len(bins_list) - 1))
            df_gwas_list.append(gwas_chr_df)
        df_input_torus = pd.concat(df_gwas_list)
        logging.info(f'prepare fastenloc gwas file all columns duration {datetime.now() - start_time}')

        # torus input file include columns variant_id、loc、zscore
        df_input_torus.to_csv(output_torus_input_file, columns=['variant_id', 'loc', 'zscore'],
                              sep=const.output_spliter, header=False, index=False)
        os.system(shell_compress_file.format(output_torus_input_file))

        #
        os.system(shell_command_torus_execute.format(f'{output_torus_input_file}.gz', output_torus_output_file))
        os.system(shell_compress_file.format(output_torus_output_file))

        logging.info(f'prepare fastenloc gwas file at: {datetime.now()}, duration： {datetime.now() - start_time}')

        return output_torus_output_file
        # clean temp file folder
        # utils.delete_dir(output_torus_input_dir)


if __name__ == '__main__':
    fastenloc_gwas_pro = FastenlocGwasProcessor()
    processor = gdp.Processor()
    _working_dir = os.path.join(processor.tool_parent_dir, 'fastenloc')
    fastenloc_gwas_pro.prepare_gwas_data(working_dir=_working_dir,
                                         gwas_preprocessed_file=processor.gwas_preprocessed_file,
                                         gwas_col_dict=processor.gwas_col_dict,
                                         var_id_col_name=gdp.Processor.VAR_ID_COL_NAME,
                                         ld_block_loci_file=processor.global_config['input']['ld_block_loci_file'])
