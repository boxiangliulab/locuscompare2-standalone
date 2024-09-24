import logging
import os
from datetime import datetime
from pathlib import Path
import sys
import pandas as pd
from common import global_data_process as gdp, constants as const, coloc_utils as utils
import ast
import concurrent
import json
import logging
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import yaml


def outputschedule(rownum, totalnum, currenttissuenum, numoftissues, rank_dir):
    calculated_schedule = int(rownum/totalnum * 80/numoftissues + 80/numoftissues * (currenttissuenum - 1))
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()


## original

# class FastenlocGwasProcessor:

#     def __init__(self):
#         logging.info('init FastenlocGwasProcessor')

#     def __loc_array_init(self, chr_str, chr_len):
#         arr = [f'loc_{chr_str}_{x}' for x in range(chr_len)]
#         return arr

#     def get_output_gwas_dir(self, output_base_dir):
#         return f'{output_base_dir}/gwas'

#     def get_output_torus_output(self, output_gwas_dir):
#         return f'{output_gwas_dir}/pip'

#     def get_output_torus_output_file(self, output_torus_output_dir):
#         return f'{output_torus_output_dir}/torus_output.pip'

#     def prepare_gwas_data(self, working_dir=None, gwas_preprocessed_file=None,
#                           gwas_col_dict=None, ld_block_loci_file=None, rank_dir=None,
#                           currenttissuenum=None, numoftissues=None, 
#                           whether_schedual=False):
#         start_time = datetime.now()
#         logging.info(f'Grouping fastenloc preprocessed gwas file at {start_time}')

#         output_base_dir = working_dir
#         output_gwas_dir = self.get_output_gwas_dir(output_base_dir)
#         output_torus_input_dir = f'{output_gwas_dir}/torus_input'
#         output_torus_input_file = f'{output_torus_input_dir}/torus_input'
#         utils.delete_file_if_exists(output_torus_input_file)
#         output_torus_output_dir = self.get_output_torus_output(output_gwas_dir)
#         output_torus_output_file = self.get_output_torus_output_file(output_torus_output_dir)

#         dap_use_col = [gwas_col_dict['beta'],
#                        gwas_col_dict['se'],
#                        gwas_col_dict['position'], gwas_col_dict['chrom'],
#                        gwas_col_dict['effect_allele'],
#                        gwas_col_dict['other_allele'], 'ref_', 'alt_']

#         shell_command_torus_execute = 'torus -d {} --load_zval -dump_pip {}'
#         shell_compress_file = 'gzip -k -f {}'

#         Path(output_torus_input_dir).mkdir(parents=True, exist_ok=True)
#         Path(output_torus_output_dir).mkdir(parents=True, exist_ok=True)
#         gwas_snp_count = 0
#         tmp = pd.read_csv(gwas_preprocessed_file, sep=const.column_spliter)
#         total_num = len(tmp)
#         del tmp
#         ix = 0
#         with pd.read_table(gwas_preprocessed_file, sep=const.column_spliter, usecols=dap_use_col,
#                            dtype={gwas_col_dict['position']: 'Int64',
#                                   gwas_col_dict['chrom']: 'category',
#                                   gwas_col_dict['effect_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
#                                   gwas_col_dict['other_allele']: pd.CategoricalDtype(const.SNP_ALLELE),
#                                   'ref_': pd.CategoricalDtype(const.SNP_ALLELE),
#                                   'alt_': pd.CategoricalDtype(const.SNP_ALLELE)},
#                            iterator=True, chunksize=500000) as reader:
#             for chunk in reader:
#                 ix = ix + len(chunk)
#                 if whether_schedual == True:
#                     outputschedule(rownum=ix,
#                                 totalnum=total_num,
#                                 currenttissuenum = currenttissuenum,
#                                 numoftissues=numoftissues,
#                                 rank_dir=rank_dir)
#                 gwas_snp_count += chunk.shape[0]
#                 chunk['zscore'] = chunk[gwas_col_dict['beta']] / chunk[gwas_col_dict['se']]
#                 chunk.drop(columns=[gwas_col_dict['beta'], gwas_col_dict['se']], inplace=True)
#                 # Adjust allele order
#                 ref_df = chunk[[gwas_col_dict['chrom'], gwas_col_dict['position'], 'alt_', 'ref_']]
#                 logging.info(f'Aligning allele order')
#                 utils.adjust_allele_order(chunk,
#                                           gwas_col_dict['effect_allele'],
#                                           gwas_col_dict['other_allele'],
#                                           gwas_col_dict['chrom'],
#                                           gwas_col_dict['position'],
#                                           ref_df,
#                                           ref_df_chrom_col_name=gwas_col_dict['chrom'],
#                                           ref_df_pos_col_name=gwas_col_dict['position'],
#                                           ref_df_alt_allele_col_name='alt_',
#                                           ref_df_ref_allele_col_name='ref_',
#                                           gz_col_name='zscore',
#                                           drop_ref_df_non_intersect_items=False)
#                 logging.info(f'Align allele order completed')
#                 del ref_df
#                 # chr1_13550_G_A_b38
#                 chunk['variant_id'] = 'chr' + chunk[gwas_col_dict['chrom']].astype(str) + '_' + \
#                                       chunk[gwas_col_dict['position']].astype(str) + '_' + \
#                                       chunk['ref_'].astype(str) + '_' + chunk['alt_'].astype(str)
#                 chunk.drop(columns=[gwas_col_dict['effect_allele'], gwas_col_dict['other_allele'], 'ref_', 'alt_'],
#                            inplace=True)
#                 logging.info(f'Adding zscore/variant_id completed, memory usage: {chunk.memory_usage(deep=True)}')
#                 # get loc range in eur_ld.hg38.bed for loc column
#                 df_ld_bed_loc = pd.read_csv(ld_block_loci_file, sep=const.column_spliter,
#                                             dtype={'start': 'Int64', 'stop': 'Int64'})
#                 df_ld_bed_loc.chr = df_ld_bed_loc.chr.apply(lambda x: int(x.replace('chr', '')))

#                 df_gwas_chr_list = list(chunk.groupby(gwas_col_dict['chrom']))
#                 df_gwas_list = []
#                 for gwas_chr in df_gwas_chr_list:
#                     chr_num = int(gwas_chr[0])
#                     gwas_chr_df = gwas_chr[1]
#                     bins_list = df_ld_bed_loc[df_ld_bed_loc['chr'] == chr_num].loc[:, 'start'] - 1
#                     gwas_chr_df['loc'] = pd.cut(x=gwas_chr_df[gwas_col_dict['position']],
#                                                 bins=list(bins_list),
#                                                 labels=self.__loc_array_init(chr_num, len(bins_list) - 1))
#                     gwas_chr_df.drop(columns=[gwas_col_dict['chrom'], gwas_col_dict['position']], inplace=True)
#                     df_gwas_list.append(gwas_chr_df)
#                 if os.path.exists(output_torus_input_file) and os.path.getsize(output_torus_input_file) > 0:
#                     mode = 'a'
#                     header = False
#                 else:
#                     mode = 'w'
#                     header = True
#                 df_input_torus = pd.concat(df_gwas_list)
#                 # torus input file include columns variant_id、loc、zscore
#                 df_input_torus.to_csv(output_torus_input_file, columns=['variant_id', 'loc', 'zscore'],
#                                       sep=const.output_spliter, mode=mode, header=header, index=False)
                
#                 del df_gwas_list, df_input_torus
#         logging.info(f'prepare fastenloc gwas file all columns duration {datetime.now() - start_time}')
#         os.system(shell_compress_file.format(output_torus_input_file))
#         # os.system(f"cp {output_torus_input_file}.gz ~/scratch/torus.input")
#         torus_cmd = shell_command_torus_execute.format(f'{output_torus_input_file}.gz', output_torus_output_file)
#         os.system(torus_cmd)
#         logging.info(f'torus cmd: {torus_cmd}')
#         os.system(shell_compress_file.format(output_torus_output_file))
#         # os.system(f"cp {output_torus_output_file}.gz ~/scratch/torus.output")
#         logging.info(f'prepare fastenloc gwas file at: {datetime.now()}, duration: {datetime.now() - start_time}')
#         return output_torus_output_file, gwas_snp_count
#         # clean temp file folder
#         # utils.delete_dir(output_torus_input_dir)




# # new version fastenloc & full gene
# class FastenlocGwasProcessor:

#     def __init__(self):
#         logging.info('init FastenlocGwasProcessor')

#     def __loc_array_init(self, chr_str, chr_len):
#         arr = [f'loc_{chr_str}_{x}' for x in range(chr_len)]
#         return arr

#     def get_output_gwas_dir(self, output_base_dir):
#         return f'{output_base_dir}/gwas'

#     def get_output_torus_output(self, output_gwas_dir):
#         return f'{output_gwas_dir}/pip'

#     def get_output_torus_output_file(self, output_torus_output_dir):
#         return f'{output_torus_output_dir}/torus_output.pip'

#     def prepare_fastenlocinput_data(self, working_dir=None, 
#                                     gwas_preprocessed_file=None,
#                                     gwas_col_dict=None, eqtl_file = None, 
#                                     eqtl_col_dict = None, rank_dir=None,
#                                     currenttissuenum=None, numoftissues=None, 
#                                     whether_schedual=False):
#         start_time = datetime.now()
#         logging.info(f'Grouping fastenloc preprocessed gwas file at {start_time}')

#         output_base_dir = working_dir
#         output_gwas_dir = self.get_output_gwas_dir(output_base_dir)
#         output_prefastenloc_dir = f'{output_gwas_dir}/fastenloc_input'
#         output_prefastenloc_file = f'{output_prefastenloc_dir}/pre_fastenloc_summarystat.tsv.gz'
#         utils.delete_file_if_exists(output_prefastenloc_file)


#         gwas_use_col = [gwas_col_dict['snp'], gwas_col_dict['beta'],
#                        gwas_col_dict['se']]
        

#         eqtl_use_col = [eqtl_col_dict['chrom'], eqtl_col_dict['gene_id'], eqtl_col_dict['snp'], eqtl_col_dict['beta'], eqtl_col_dict['se']]
        
#         eqtl_df = pd.read_csv(eqtl_file, sep=const.column_spliter, usecols=eqtl_use_col)

#         Path(output_prefastenloc_dir).mkdir(parents=True, exist_ok=True)


#         gwas_preprocessed_df = pd.read_table(gwas_preprocessed_file, sep=const.column_spliter, usecols=gwas_use_col)


#         eqtl_df.index = eqtl_df[eqtl_col_dict['se']]
#         gwas_preprocessed_df.index = gwas_preprocessed_df[gwas_col_dict['snp']]



#         # # 获取共同的 SNP
#         common_snp = list(set(gwas_preprocessed_df.index) & set(eqtl_df.index))
#         # # 筛选共同的 SNP
#         eqtl_df = eqtl_df.loc[common_snp]
#         gwas_preprocessed_df = gwas_preprocessed_df.loc[common_snp]
#         # # 合并 eqtl 和 gwas 数据

#         merged_df = pd.merge(eqtl_df[[eqtl_col_dict['gene_id'], eqtl_col_dict['snp'], eqtl_col_dict['beta'], eqtl_col_dict['se']]], 
#                              gwas_preprocessed_df[[gwas_col_dict['beta'], gwas_col_dict['se']]], 
#                              left_index=True, 
#                              right_index=True, 
#                              how='left')

#         # 保存合并的数据到指定路径
#         merged_df.to_csv(output_prefastenloc_file, sep=const.output_spliter, header=False, index=False, compression='gzip')

        
#         if whether_schedual == True:
#             outputschedule(currenttissuenum = currenttissuenum,
#                         numoftissues=numoftissues,
#                         rank_dir=rank_dir)

#         # os.system(f"cp {output_torus_output_file}.gz ~/scratch/torus.output")
#         logging.info(f'prepare fastenloc gwas file at: {datetime.now()}, duration: {datetime.now() - start_time}')
#         return output_prefastenloc_file



# new version fastenloc & gwas loci
class FastenlocGwasProcessor:
    def __init__(self):
        logging.info('init FastenlocGwasProcessor')

    def __loc_array_init(self, chr_str, chr_len):
        arr = [f'loc_{chr_str}_{x}' for x in range(chr_len)]
        return arr

    def get_output_preprocess_dir(self, output_base_dir):
        return f'{output_base_dir}/fastenloc_input'


    def __convert_positions_str_to_list(self, positions_str):

        if isinstance(positions_str, str):
            return json.loads(positions_str)
        elif isinstance(positions_str, list):
            return positions_str
        else:
            return []

    def __get_cluster_significant_snps_dict(self, cluster_df):

        cluster_snps_dict = {}
        for _, row in cluster_df.iterrows():
            cluster_snps_dict[row.loc['range_lead']] = self.__convert_positions_str_to_list(row.loc['positions'])
        return cluster_snps_dict
    
    def prepare_fastenlocinput_data(self, 
                                    working_dir=None,
                                    var_id_col_name=None,
                                    gwas_cluster_output_dir=None,
                                    gwas_cluster_summary=None,
                                    gwas_col_dict=None,
                                    eqtl_output_report=None,
                                    eqtl_output_dir=None,
                                    eqtl_col_dict=None,
                                    parallel=False,
                                    parallel_worker_num=10, 
                                    rank_dir = None,
                                    currenttissuenum = None, 
                                    numoftissues = None, 
                                    whether_schedual = False):
        


        start_time = datetime.now()
        Path(working_dir).mkdir(parents=True, exist_ok=True)
        logging.info(f'Grouping fastenloc preprocessed gwas file at {start_time}')

        eqtl_summary_df = pd.read_csv(eqtl_output_report, sep=const.column_spliter,
                                      dtype={eqtl_col_dict['chrom']: 'category'})
        gwas_summary_df = pd.read_csv(gwas_cluster_summary, sep=const.column_spliter,
                                      dtype={gwas_col_dict['chrom']: 'category'})
        gwas_cluster_snps_dict = self.__get_cluster_significant_snps_dict(gwas_summary_df)
        del gwas_summary_df
        # Put gwas range files in a list
        gwas_range_files = {}
        for gwas_range_file in os.listdir(gwas_cluster_output_dir):
            part_list = utils.get_file_name(gwas_range_file).split('-')
            if len(part_list) < 2 or 'chr' not in part_list[1]:
                continue
            chrom = part_list[1].strip('chr')
            range_files = gwas_range_files.get(chrom, [])
            range_files.append(os.path.join(gwas_cluster_output_dir, gwas_range_file))
            gwas_range_files[chrom] = range_files

        # Loop to process all eQTL trait file
        gwas_chroms = gwas_range_files.keys()

        total_len = len(eqtl_summary_df)

        output_base_dir = working_dir
        output_preprocess_dir = self.get_output_preprocess_dir(output_base_dir)
        output_prefastenloc_file = f'{output_preprocess_dir}/pre_fastenloc_summarystat.tsv.gz'
        Path(output_preprocess_dir).mkdir(parents=True, exist_ok=True)
        utils.delete_file_if_exists(output_prefastenloc_file)
        # if parallel:
        #     with ThreadPoolExecutor(max_workers=parallel_worker_num) as executor:
        #         futures = []
        #         for ix, row in eqtl_summary_df.iterrows():
        #             if whether_schedual:
        #                 outputschedule(rownum=ix,
        #                     totalnum=total_len,
        #                     currenttissuenum = currenttissuenum,
        #                     numoftissues=numoftissues,
        #                     rank_dir=rank_dir)
                        
        #             chrom = str(row.loc['chrom'])
        #             if chrom not in gwas_chroms:
        #                 continue
        #             eqtl_gene_file = os.path.join(eqtl_output_dir, chrom, row.loc['gene_file'])
        #             gene_id = utils.get_file_name(eqtl_gene_file)
        #             for gwas_range_file in gwas_range_files[chrom]:
        #                 eqtl_significant_positions = ast.literal_eval(row.loc['positions'])
        #                 range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
        #                 if len(set(gwas_cluster_snps_dict[range_lead_snp]) & set(eqtl_significant_positions)) == 0:
        #                     continue
        #                 futures.append(executor.submit(self.process_gene, working_dir,
        #                                                gwas_range_file, gwas_col_dict, row,
        #                                                eqtl_gene_file,
        #                                                var_id_col_name, gene_id, eqtl_col_dict, output_prefastenloc_file))

        #         for future in concurrent.futures.as_completed(futures):
        #             try:
        #                 data = future.result()
        #             except Exception as exc:
        #                 logging.error("".join(traceback.TracebackException.from_exception(exc).format()))

        # else:
        for ix, row in eqtl_summary_df.iterrows():
            if whether_schedual:
                outputschedule(rownum=ix,
                    totalnum=total_len,
                    currenttissuenum = currenttissuenum,
                    numoftissues=numoftissues,
                    rank_dir=rank_dir)
                    
            chrom = str(row.loc['chrom'])
            if chrom not in gwas_chroms:
                continue
            eqtl_gene_file = os.path.join(eqtl_output_dir, chrom, row.loc['gene_file'])
            gene_id = utils.get_file_name(eqtl_gene_file)
            for gwas_range_file in gwas_range_files[chrom]:
                eqtl_significant_positions = ast.literal_eval(row.loc['positions'])
                range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
                if len(set(gwas_cluster_snps_dict[range_lead_snp]) & set(eqtl_significant_positions)) == 0:
                    continue
                self.process_gene(working_dir, gwas_range_file,
                                    gwas_col_dict, row, eqtl_gene_file,
                                    var_id_col_name, gene_id, eqtl_col_dict, output_prefastenloc_file)



        logging.info(f'prepare fastenloc gwas file at: {datetime.now()}, duration: {datetime.now() - start_time}')
        return output_prefastenloc_file



    def process_gene(self, working_dir, gwas_range_file, gwas_col_dict, row, eqtl_gene_file,
                     var_id_col_name, gene_id, eqtl_col_dict, output_prefastenloc_file):

        range_lead_snp = utils.get_file_name(gwas_range_file).split('-')[0]
        candidate_gwas_df = pd.read_table(gwas_range_file, sep=const.column_spliter)
        if len(candidate_gwas_df) <= 1:
            return
        eqtl_trait_df = pd.read_table(eqtl_gene_file, sep=const.column_spliter)

        eqtl_trait_df.drop(
            index=eqtl_trait_df[~eqtl_trait_df[var_id_col_name].isin(candidate_gwas_df[var_id_col_name])].index,
            inplace=True)
        if len(eqtl_trait_df) <= 1:
            return
        utils.drop_non_intersect_rows(eqtl_trait_df, var_id_col_name, candidate_gwas_df, var_id_col_name)
        if len(candidate_gwas_df) <= 1:
            return


        # candidate_gwas_df.reset_index(drop=True, inplace=True)
        # eqtl_trait_df.reset_index(drop=True, inplace=True)

        candidate_gwas_df.index = candidate_gwas_df[var_id_col_name]
        eqtl_trait_df.index = eqtl_trait_df[var_id_col_name]
        eqtl_trait_df['ix'] = f'{gene_id}_{range_lead_snp}'

        merged_df = pd.merge(eqtl_trait_df[['ix', eqtl_col_dict['snp'], eqtl_col_dict['beta'], eqtl_col_dict['se']]], 
                             candidate_gwas_df[[gwas_col_dict['beta'], gwas_col_dict['se']]], 
                             left_index=True, 
                             right_index=True, 
                             how='left')


        if os.path.exists(output_prefastenloc_file) and os.path.getsize(output_prefastenloc_file) > 0:
            mode = 'a'
            header = False
        else:
            mode = 'w'
            header = False

        # torus input file include columns variant_id、loc、zscore
        merged_df.to_csv(output_prefastenloc_file, 
                                sep=const.output_spliter, mode=mode, header=header, index=False)


if __name__ == '__main__':
    fastenloc_gwas_pro = FastenlocGwasProcessor()
    processor = gdp.Processor()
    _working_dir = os.path.join(processor.tool_parent_dir, 'fastenloc')
    fastenloc_gwas_pro.prepare_gwas_data(working_dir=_working_dir,
                                         gwas_preprocessed_file=processor.gwas_preprocessed_file,
                                         gwas_col_dict=processor.gwas_col_dict,
                                         ld_block_loci_file=processor.global_config['input']['ld_block_loci_file'])
