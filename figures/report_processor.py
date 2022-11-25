import pandas as pd
import json
import os
from pathlib import Path
from common import coloc_utils as utils, constants as const
from figures import plot_processsor as plotp
from datetime import datetime


class ReportProcessor:
    def __init__(self, global_processor):
        print('init ReportDataProcessor')
        plot = plotp.Plot(global_processor)
        self.plot = plot
        self.coloc_gwas_merge_from_preprocessed_file = True
        self.gwas_preprocessed_file = global_processor.gwas_preprocessed_file
        self.tool_parent_dir = global_processor.tool_parent_dir
        self.gwas_preprocessed_file_pd = None
        self.report_col = ['gwas_path', 'eqtl_path', 'gene_id', 'rsid', 'chrom']
        self.report_head_ten_for_tool = {'origin': {}, 'filter_rsid': {}, 'filter_gene': {}, 'manht': {},
                                         'profile': {'tissue': global_processor.eqtl_tissue,
                                                     "trait": global_processor.gwas_trait}}
        self.output_dir_path = f'{self.tool_parent_dir}/figures/'

        self.input_read_gwas_columns = plot.input_read_gwas_columns

        Path(self.output_dir_path).mkdir(parents=True, exist_ok=True)

    # 存储每个tool的report数据
    def init_reprot_data(self, report_file_path=None, report_df=None, tool_name=None, lens=10, report_col=None):
        print(report_file_path)
        if not utils.file_exists(report_file_path):
            return
        if report_col is None:
            report_col = self.report_col
        if report_df is None:

            if tool_name == 'coloc':
                report_df = pd.read_csv(report_file_path, sep=const.column_spliter).groupby(
                    'gene_id').max().reset_index().sort_values(by=['overall_H4', 'SNP.PP.H4'], ascending=False,
                                                               inplace=False)[report_col]
                print(report_df.columns)
            else:
                report_df = pd.read_csv(report_file_path, sep=const.column_spliter, usecols=report_col)
        self.report_head_ten_for_tool['origin'][tool_name] = report_df.head(lens)

    # 处理每组coloc plot的三个图的文件数据，生成jsonp文件
    def prepare_jsonp_data_(self):
        # 处理所有工具数据对应的gwas和eqtl josnp文件，重复文件的不处理
        print('origin:', self.report_head_ten_for_tool['origin'])
        for col_name in self.report_head_ten_for_tool['origin']:
            print(col_name)
            for v, col in self.report_head_ten_for_tool['origin'][col_name].iterrows():
                # 相同的snp的gene_id文件需要一起merge，然后生成文件，因为最终绘图，相同snp的gene_id需要一起展示
                # 如果snp 对应的 gene_id数组里面已经有该基因了 则不累加
                self.filter_list_map(self.report_head_ten_for_tool['filter_rsid'], col, 'rsid', 'gene_id', col_name, v)
                # self.filter_list_map(self.report_head_ten_for_tool['filter_gene'], col, 'gene_id', 'rsid')

            self.report_head_ten_for_tool['origin'][col_name] = self.report_head_ten_for_tool['origin'][
                col_name].to_dict(orient='records')

        if self.coloc_gwas_merge_from_preprocessed_file:
            print('merge preprocessed file')
            for rsid_name in self.report_head_ten_for_tool['filter_rsid'].keys():
                rsid_to_gene_obj = self.report_head_ten_for_tool['filter_rsid'][rsid_name]
                self.plot.coloc_range_gwas_preprocessed_merge_eqtl(self.gwas_preprocessed_file_pd, rsid_to_gene_obj,
                                                                   self.output_dir_path,
                                                                   rsid_name)
        else:
            print('merge clustered file')
            for rsid_name in self.report_head_ten_for_tool['filter_rsid'].keys():
                rsid_to_gene_obj = self.report_head_ten_for_tool['filter_rsid'][rsid_name]
                print(rsid_to_gene_obj)
                merge_pd = self.plot.prepare_gwas_map_gene_list_df(rsid_to_gene_obj)
                if merge_pd is None:
                    continue
                self.plot.coloc_range_gwas_merge_eqtl(merge_pd, rsid_to_gene_obj, self.output_dir_path, rsid_name)

        # 将所有工具总体数据生成jsonp文件用于report展示
        output_file = f'{self.output_dir_path}/data/tools_report_data.json'
        with open(output_file, 'w', encoding='utf-8') as output_file:
            output_file.write(
                'callbackGetToolsReportData(' + json.dumps(self.report_head_ten_for_tool, indent=2) + ')')

    # 将report 总数据重新组合成 rsid to gene mapping
    def filter_list_map(self, filter_list, col, type_str, value_str, tool, index):
        key_value = col[type_str]
        # predixcan report 没有rsid 也没有var_id 暂时使用 tool_rank_index 代替
        if pd.isnull(key_value):
            key_value = f'{tool}_rank_{index + 1}'

        item_value = col[value_str]
        if filter_list.get(key_value) is None:
            filter_list[key_value] = {'path': col.gwas_path, 'list': {}, 'chrom': col.chrom}

        is_exist = False
        for item_id in filter_list[key_value]['list'].keys():
            if item_id == item_value:
                is_exist = True
                break
        if is_exist == False:
            filter_list[key_value]['list'][item_value] = {"path": col.eqtl_path}
        return filter_list

    # 处理manhattan plot数据，生成jsonp文件
    def show_manhattan_plot(self):
        # gwas_preprocessed_file only read once in report precossor
        self.gwas_preprocessed_file_pd = pd.read_table(self.gwas_preprocessed_file, sep=const.column_spliter, header=0,
                                                       usecols=self.input_read_gwas_columns)[
            self.input_read_gwas_columns]
        self.report_head_ten_for_tool['manht']['path'] = self.plot.show_manhattan_plot(self.output_dir_path,
                                                                                       self.gwas_preprocessed_file_pd)

    # 处理manhattan 数据
    # 处理coloc 数据
    # 处理js文件
    def to_show_html_page(self):
        start_time = datetime.now()
        if len(self.report_head_ten_for_tool['origin'].keys()) == 0:
            return
        self.show_manhattan_plot()
        self.prepare_jsonp_data_()
        # 处理html，js文件模版，生成对应的前端文件夹
        self.plot.report_html_handler(self.output_dir_path)
        print(f'report complete at: {datetime.now()},duration: {datetime.now() - start_time}')

    def find_output_report_file(self, tool_name):
        file_path = f'{self.tool_parent_dir}/{tool_name}/analyzed'
        for output_file_name in os.listdir(file_path):
            if output_file_name.startswith(f'{tool_name}_output'):
                return f'{file_path}/{output_file_name}'

# if __name__ == '__main__':
#     gwas_trait = 'T2D'
#     eqtl_list = ['Adipose_Visceral_Omentum', 'Colon_Sigmoid', 'Colon_Transverse', 'Kidney_Cortex', 'Liver',
#                  'Muscle_Skeletal', 'Pancreas', 'Small_Intestine_Terminal_Ileum', 'Stomach']
#     for tissue in eqtl_list:
#         tool_parent_dir = f'/Volumes/HD/biodata/colocalization-tools/processed/EUR_T2D/{gwas_trait}-{tissue}'
#         if utils.file_exists(tool_parent_dir):
#             rp = ReportProcessor(
#                 tool_parent_dir=tool_parent_dir,
#                 gwas_preprocessed_file='/Volumes/HD/biodata/colocalization-tools/preprocessed/gwas/T2D/preprocessed.tsv',
#                 eqtl_tissue=tissue, gwas_trait='T2D')
#             rp.init_reprot_data(report_file_path=rp.find_output_report_file('jlim'), tool_name='jlim', lens=10)
#             rp.init_reprot_data(report_file_path=rp.find_output_report_file('fastenloc'), tool_name='fastenloc',
#                                 lens=10)
#             rp.init_reprot_data(report_file_path=rp.find_output_report_file('coloc'), tool_name='coloc', lens=10)
#             rp.init_reprot_data(report_file_path=rp.find_output_report_file('smr'), tool_name='smr', lens=10)
#             rp.init_reprot_data(report_file_path=rp.find_output_report_file('predixcan'), tool_name='predixcan',
#                                 lens=10)
#             rp.to_show_html_page()
