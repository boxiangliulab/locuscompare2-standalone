"""
Process GWAS Catalog harmonised file in http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
Columns in this file:
hm_variant_id	hm_rsid hm_chrom	hm_pos	hm_other_allele	hm_effect_allele	hm_beta	hm_odds_ratio
hm_ci_lower	hm_ci_upper	hm_effect_allele_frequency	hm_code	chromosome	base_pair_location	effect_allele
other_allele	n	effect_allele_frequency	beta	standard_error	p_value	odds_ratio	ci_lower	ci_upper	variant_id
"""

import datetime
import pandas as pd
from common import constants as const, coloc_utils as utils
from pathlib import Path
import json
import os
import shutil
import re
from figures import sql_query_snp_r2 as sqr
import logging


class Plot:
    def __init__(self, processor):
        logging.info('init Plot')
        self.gwas_col_dict = processor.gwas_col_dict
        self.eqtl_col_dict = processor.eqtl_col_dict
        self.var_id = processor.VAR_ID_COL_NAME
        self.population = processor.global_config.get('population', 'EUR').upper()
        # 读取输入文件所需列名
        self.input_read_gwas_columns = [self.gwas_col_dict["snp"], self.gwas_col_dict["chrom"],
                                        self.gwas_col_dict["position"],
                                        self.gwas_col_dict["pvalue"], self.var_id]
        self.input_read_eqtl_columns = [self.eqtl_col_dict["snp"], self.eqtl_col_dict["chrom"],
                                        self.eqtl_col_dict["position"],
                                        self.eqtl_col_dict["pvalue"], self.var_id]
        # josnp 输出文件列名
        self.output_read_columns_mode = ['snp', 'chrom', 'position', 'pvalue', 'var_id', 'category']
        self.output_read_columns = {
            'snp': 'snp',
            'chrom': 'chrom',
            'position': 'position',
            'pvalue': 'pvalue',
            'var_id': 'var_id',
            'category': 'category'
        }

    def filter_data_frame_by_p_value(self, input_df, p_value, pval_max_count, inplace=True):
        p_value_col_name = self.output_read_columns['pvalue']
        bind_data_df = input_df  # 存储有效数据
        bins_list = [p_value]

        # drop rows with p-value > p_value_threshold until data less  pval_max_count
        while len(bind_data_df) > pval_max_count:
            p_value = p_value * 0.1
            bind_data_df = utils.filter_data_frame_by_p_value(bind_data_df, p_value, p_value_col_name, inplace)
            bins_list.insert(0, p_value)

        # 根据上面截取数据的pvalue，随机获取大于pvalue的1w条数据作为画图数据
        # 这1w条数据根据数组 bins_list 分段随机获取等比例的数据，这样在manhattan图上下面的点就不会出现断层
        # 例如 bins_list= [0.0005, 0.005, 0.05, 0.5] ,则在0～0.0005、0.0005~0.005、0.005~0.05、0.05~0.5四个区间分别随机拿2.5k条组成1w条用于填充图下方空白区域
        logging.debug(f'sample data params :  {p_value} {bins_list}')
        unbind_data_df = input_df.drop(
            input_df[(input_df[p_value_col_name] < p_value) | (input_df[p_value_col_name] == 0.0)].index,
            inplace=inplace)
        p_value_group = pd.cut(unbind_data_df[p_value_col_name], bins=bins_list)
        unbind_data_df = unbind_data_df.groupby(p_value_group).sample(n=round(10000 / (len(bins_list) - 1)),
                                                                      random_state=1)
        logging.debug(f'after parse bind data lines: {len(bind_data_df)} , unbind data lines: {len(unbind_data_df)} ')

        return pd.concat([bind_data_df, unbind_data_df])

    # 筛选input file中pvalue字段值小等于manhattan_min_pval的数据
    def __parse_tsv_data_to_json(self, pval, pval_max_count, output_file, gwas_df, output_dir):

        df_len = len(gwas_df)

        # utils.clean_data(gwas_df, dup_consider_subset="variant_id")
        if df_len < pval_max_count:
            pval_filter_gwas_df = gwas_df.drop(gwas_df[gwas_df[self.output_read_columns['pvalue']] == 0.0].index,
                                               inplace=False)
        else:
            pval_filter_gwas_df = self.filter_data_frame_by_p_value(gwas_df, pval, pval_max_count, inplace=False)
        logging.info(f'parse finish, current file lines: {len(pval_filter_gwas_df)}, output_file:{output_file}')

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir + '/data/').mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', encoding='utf-8') as output_file:
            output_file.write(
                'callbackGetPlotData(' + json.dumps(pval_filter_gwas_df.to_dict(orient='records'), indent=2) + ')')

    def __file_path_parse(self, file_pd, output_dir, title_name, type_str):
        out_file_path = output_dir + 'data/' + type_str + '_' + title_name + '.json'
        # 重复文件的不处理
        if not utils.file_exists(out_file_path):
            self.__parse_tsv_data_to_json(const.manhattan_min_pval, const.manhattan_max_plot, out_file_path, file_pd,
                                          output_dir)
        return 'data/' + type_str + '_' + title_name + '.json'

    def report_html_handler(self, output_dir):
        # 复制项目static目录文件到所需展示数据目录下
        current_dir = os.path.dirname(Path(__file__).resolve())
        static_files_path = os.path.join(current_dir, 'static')
        output_static_files_path = os.path.join(output_dir, 'static')
        self.__copy_static_files(static_files_path, output_static_files_path)

        # 复制html文件所需展示数据目录下
        shutil.copy(os.path.join(current_dir, 'index.html'),
                    output_dir + 'report.html')
        shutil.copy(os.path.join(static_files_path, 'mode.js'),
                    output_static_files_path + '/report.js')

    def __copy_static_files(self, source_path, target_path):
        if not os.path.exists(target_path):
            os.makedirs(target_path)

        if os.path.exists(source_path):
            # root 所指的是当前正在遍历的这个文件夹的本身的地址
            # dirs 是一个 list，内容是该文件夹中所有的目录的名字(不包括子目录)
            # files 同样是 list, 内容是该文件夹中所有的文件(不包括子目录)
            for root, dirs, files in os.walk(source_path):
                for file in files:
                    src_file = os.path.join(root, file)
                    shutil.copy(src_file, target_path)

    def coloc_range_gwas_preprocessed_merge_eqtl(self, gwas_preprocessed_pd, path_obj, output_dir, rsid):
        path_obj['path'] = None
        for gene_id in path_obj['list'].keys():
            eqtl_file_path = path_obj['list'][gene_id]['path']
            if eqtl_file_path is None or pd.isnull(eqtl_file_path) or not utils.file_exists(eqtl_file_path):
                logging.info(f'file no exist {eqtl_file_path}')
                path_obj['list'][gene_id]['path'] = None
                path_obj['list'][gene_id]['gwas_path'] = None
            else:
                eqtl_file_path_pd = pd.read_table(eqtl_file_path, sep=const.column_spliter, header=0,
                                                  usecols=self.input_read_eqtl_columns)
                eqtl_file_path_pd.rename(
                    columns={self.eqtl_col_dict["snp"]: 'eqtl_snp',
                             self.eqtl_col_dict["position"]: 'eqtl_position',
                             self.eqtl_col_dict["chrom"]: 'eqtl_chrom',
                             self.eqtl_col_dict["pvalue"]: 'eqtl_pvalue',
                             },
                    inplace=True)
                merge_pd = pd.merge(left=gwas_preprocessed_pd,
                                    right=eqtl_file_path_pd,
                                    how='inner', on=self.var_id)
                r2_pd = sqr.retrieve_ld(chr=path_obj['chrom'], rsid=rsid, population=self.population,
                                        rsid_col_name=self.gwas_col_dict["snp"])
                merge_pd = pd.merge(left=merge_pd, right=r2_pd, how='left', on=self.gwas_col_dict["snp"])
                merge_pd['r2'].fillna(0, inplace=True)
                merge_pd['category'] = pd.cut(x=merge_pd['r2'], bins=[0, 0.2, 0.4, 0.6, 0.8, 1, 2],
                                              labels=['r2_one', 'r2_two', 'r2_thr', 'r2_four', 'r2_five',
                                                      'target_pvalue'])
                merge_pd['category'].fillna('r2_one', inplace=True)

                if len(merge_pd) == 0:
                    path_obj['list'][gene_id]['path'] = None
                    path_obj['path'] = None
                else:
                    # crear gwas_rsid_gene.json file
                    gwas_merge_pd = merge_pd[['eqtl_snp', self.gwas_col_dict["chrom"],
                                        self.gwas_col_dict["position"],
                                        self.gwas_col_dict["pvalue"], self.var_id] + ['category']]
                    gwas_merge_pd = gwas_merge_pd.rename(
                        columns={'eqtl_snp': self.output_read_columns['snp'],
                                 self.gwas_col_dict["position"]: self.output_read_columns['position'],
                                 self.gwas_col_dict["chrom"]: self.output_read_columns['chrom'],
                                 self.gwas_col_dict["pvalue"]: self.output_read_columns['pvalue'],
                                 self.var_id: self.output_read_columns['var_id'],
                                 })
                    path_obj['list'][gene_id]['gwas_path'] = self.__file_path_parse(gwas_merge_pd, output_dir,
                                                                                    re.sub(r'\W', "_", rsid) + '_' + gene_id,
                                                                                    'gwas')
                    # crear eqtl_rsid_gene.json file
                    eqtl_merge_pd = merge_pd[
                        ['eqtl_snp', 'eqtl_position', 'eqtl_chrom', 'eqtl_pvalue', self.var_id, 'category']]
                    eqtl_merge_pd = eqtl_merge_pd.rename(columns={'eqtl_snp': self.output_read_columns['snp'],
                                                                  'eqtl_position': self.output_read_columns['position'],
                                                                  'eqtl_chrom': self.output_read_columns['chrom'],
                                                                  'eqtl_pvalue': self.output_read_columns['pvalue'],
                                                                  self.var_id: self.output_read_columns['var_id']
                                                                  })
                    path_obj['list'][gene_id]['path'] = self.__file_path_parse(eqtl_merge_pd, output_dir,
                                                                               re.sub(r'\W', "_", rsid) + '_' + gene_id, 'eqtl')
                    # crear pvalue_rsid_gene.json file
                    merge_pvalue_pd = merge_pd[
                        [self.var_id, self.gwas_col_dict['pvalue'], 'eqtl_pvalue', self.gwas_col_dict["chrom"],
                         'eqtl_snp', 'category']]
                    merge_pvalue_pd = merge_pvalue_pd.rename(
                        columns={'eqtl_snp': self.output_read_columns['snp'],
                                 self.gwas_col_dict['pvalue']: self.output_read_columns['position'],
                                 self.gwas_col_dict["chrom"]: self.output_read_columns['chrom'],
                                 'eqtl_pvalue': self.output_read_columns['pvalue'],
                                 self.var_id: self.output_read_columns['var_id']
                                 })
                    merge_pvalue_pd.sort_values('position', inplace=True)
                    path_obj['list'][gene_id]['pvalue_path'] = self.__file_path_parse(merge_pvalue_pd,
                                                                                      output_dir,
                                                                                      re.sub(r'\W', "_", rsid) + '_' + gene_id, 'pvalue')

    def coloc_range_gwas_merge_eqtl(self, merge_pds, path_obj, output_dir, rsid):
        gwas_file_path = path_obj['path']
        if gwas_file_path is None:
            path_obj['path'] = None
        else:
            gwas_file_path_pd = pd.read_table(gwas_file_path, sep=const.column_spliter, header=0,
                                              usecols=self.input_read_gwas_columns)
            gwas_merge_pd = pd.merge(left=merge_pds, right=gwas_file_path_pd, how='left', on=self.var_id)
            gwas_merge_pd.loc[:, 'category'] = ''
            gwas_merge_pd.rename(columns={self.gwas_col_dict["snp"]: self.output_read_columns['snp'],
                                          self.gwas_col_dict["position"]: self.output_read_columns['position'],
                                          self.gwas_col_dict["chrom"]: self.output_read_columns['chrom'],
                                          self.gwas_col_dict["pvalue"]: self.output_read_columns['pvalue'],
                                          self.var_id: self.output_read_columns['var_id'],
                                          },
                                 inplace=True)
            gwas_merge_pd.loc[gwas_merge_pd[self.output_read_columns['snp']] == rsid, ('category')] = 'target_pvalue'
            path_obj['path'] = self.__file_path_parse(gwas_merge_pd, output_dir, rsid, 'gwas')

        for gene_id in path_obj['list'].keys():
            eqtl_file_path = path_obj['list'][gene_id]['path']
            if eqtl_file_path is None:
                path_obj['list'][gene_id]['path'] = None
            else:
                eqtl_file_path_pd = pd.read_table(eqtl_file_path, sep=const.column_spliter, header=0,
                                                  usecols=self.input_read_eqtl_columns)
                eqtl_merge_pd = pd.merge(left=merge_pds, right=eqtl_file_path_pd, how='left', on=self.var_id)
                eqtl_merge_pd.loc[:, 'category'] = ''
                eqtl_merge_pd.rename(columns={self.eqtl_col_dict["snp"]: self.output_read_columns['snp'],
                                              self.eqtl_col_dict["position"]: self.output_read_columns['position'],
                                              self.eqtl_col_dict["chrom"]: self.output_read_columns['chrom'],
                                              self.eqtl_col_dict["pvalue"]: self.output_read_columns['pvalue'],
                                              self.var_id: self.output_read_columns['var_id'],
                                              },
                                     inplace=True)
                eqtl_merge_pd.loc[
                    eqtl_merge_pd[self.output_read_columns['snp']] == rsid, ('category')] = 'target_pvalue'
                path_obj['list'][gene_id]['path'] = self.__file_path_parse(eqtl_merge_pd, output_dir,
                                                                           rsid + '_' + gene_id, 'eqtl')
        return path_obj

    def prepare_gwas_map_gene_list_df(self, path_obj):
        if not utils.file_exists(path_obj['path']):
            path_obj['path'] = None
            logging.warning(f'gwas: {path_obj["path"]} not exist')
            return None
        else:
            gwas_file_path_pd = pd.read_table(path_obj['path'], sep=const.column_spliter, header=0,
                                              usecols=[self.var_id])
            merge_pd = gwas_file_path_pd
            for gene_id in path_obj['list'].keys():
                if not utils.file_exists(path_obj['list'][gene_id]['path']):
                    path_obj['list'][gene_id]['path'] = None
                    logging.warning(f'eqtl: {path_obj["path"]} not exist')
                else:
                    eqtl_file_path_pd = pd.read_table(path_obj['list'][gene_id]['path'], sep=const.column_spliter,
                                                      header=0,
                                                      usecols=[self.var_id])
                    merge_pd = pd.merge(left=merge_pd, right=eqtl_file_path_pd, how='inner', on=self.var_id)

            return merge_pd

    def show_manhattan_plot(self, output_dir, gwas_preprocessed_file_pd):
        start_time = datetime.datetime.now()
        logging.info(f"show plot start at {start_time}")
        file_df = gwas_preprocessed_file_pd.rename(columns={self.gwas_col_dict["snp"]: self.output_read_columns['snp'],
                                                            self.gwas_col_dict["position"]: self.output_read_columns[
                                                                'position'],
                                                            self.gwas_col_dict["chrom"]: self.output_read_columns[
                                                                'chrom'],
                                                            self.gwas_col_dict["pvalue"]: self.output_read_columns[
                                                                'pvalue'],
                                                            self.var_id: self.output_read_columns['var_id'],
                                                            },
                                                   inplace=False)
        file_df.loc[:, 'category'] = ''
        file_df.loc[(file_df[self.output_read_columns['chrom']] % 2 == 0), 'category'] = 'double_chr'
        file_df.loc[(file_df[self.output_read_columns['chrom']] % 2 > 0), 'category'] = 'single_chr'
        out_file_path = self.__file_path_parse(file_df, output_dir, 'mant', 'data')

        end_time = datetime.datetime.now()
        duration = end_time - start_time
        logging.info(f"finish process in: {duration}")
        return out_file_path
