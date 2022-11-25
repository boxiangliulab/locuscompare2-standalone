from pathlib import Path
import os
import pandas as pd
from common import global_data_process as gdp, constants as const, coloc_utils as utils, config as cf
from figures import sql_query_snp_r2 as sqr
from ranking import birra as br
import logging
from datetime import datetime
import shutil
import json

'''
    predixcan report no snp column
    coloc、jlim、smr、predixcan,ecaviar report have snp column
    all tools target snp is gwas(pvalue)+eqtl(pvalue) is min
'''

output_read_columns_mode = {
    "snp": "snp",
    "chrom": "chrom",
    "position": "position",
    "pvalue": "pvalue",
}


def report_data_process(report_list):
    start_time = datetime.now()
    logging.info('start report data process ')
    # create traits file
    create_trait_file(report_list)
    logging.info(f'report complete at: {datetime.now()},duration: {datetime.now() - start_time}')


def create_trait_file(report_list=None, ):
    logging.info('create_trait_file ')
    report_list_pd = pd.DataFrame(list(report_list),
                                  columns=['trait', 'study1', 'tool_name', 'tissue', 'report_path', 'cfg_pro'])
    df_trait_list = list(report_list_pd.groupby('trait'))
    output_base_dir = None
    for trait_item in df_trait_list:
        trait_name = trait_item[0]
        trait_df = trait_item[1]

        # every trait init config param
        processor = trait_df['cfg_pro'].iloc[0]
        report_col = ['gwas_path', 'eqtl_path', 'gene_id', 'rsid', 'chrom']
        gwas_preprocessed_file = processor.gwas_preprocessed_file
        gwas_col_dict = processor.gwas_col_dict
        eqtl_col_dict = processor.eqtl_col_dict
        population = processor.global_config.get('population', 'EUR').upper()
        output_base_dir = os.path.join(processor.report_path, 'figures', 'data')
        Path(output_base_dir).mkdir(parents=True, exist_ok=True)
        trait_path = os.path.join(output_base_dir, trait_name)
        Path(trait_path).mkdir(parents=True, exist_ok=True)

        report_html_handler(os.path.join(processor.report_path, 'figures'))
        # get ranking info
        rpt = {}
        df_genes_info_path_list = []
        for tool_index, tool_row in trait_df.iterrows():
            rpt[f'{tool_row["tool_name"]}'] = tool_row['report_path']
        ranking_output_file_path = br.run_ranking(
            f'{trait_path}/emsemble_ranking_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv', rpt)
        report_gene_ranking = const.default_report_gene_ranking
        ranking_output_file_pd = pd.read_csv(ranking_output_file_path, sep=const.column_spliter,
                                             usecols=['gene_id', 'birra_ranking']).head(report_gene_ranking)

        # gene ranking mapping every tool info
        for tool_index, tool_row in trait_df.iterrows():
            tool_name = tool_row['tool_name']
            report_file_path = tool_row['report_path']
            report_df = pd.read_csv(report_file_path, sep=const.column_spliter, usecols=report_col)
            gene_report_merge_pd = pd.merge(left=ranking_output_file_pd, right=report_df, how='inner', on='gene_id')
            gene_report_merge_pd['tool_name'] = tool_name
            # get gene target snp first
            gene_report_merge_pd.drop_duplicates('gene_id', keep='first', inplace=True)
            df_genes_info_path_list.append(gene_report_merge_pd)
        df_genes_info_path_pd = pd.concat(df_genes_info_path_list)

        # create genes file
        create_gene_file(df_genes_info_path_pd, gwas_preprocessed_file, gwas_col_dict, trait_path,
                         eqtl_col_dict, population)

    # create all need to show report file
    with open(os.path.join(output_base_dir, f'reports_path.json'), 'w', encoding='utf-8') as output_file:
        output_file.write(
            'callbackGetToolsReportData(' + json.dumps(
                report_list_pd[['trait', 'tool_name', 'study1', 'tissue', 'report_path']].to_dict(
                    orient='records'), indent=2) + ')')


def create_gene_file(df_genes_info_path_pd=None, gwas_preprocessed_file=None, gwas_col_dict=None, trait_path=None,
                     eqtl_col_dict=None, population=None):
    logging.info('create_gene_file')
    var_id = gdp.Processor.VAR_ID_COL_NAME
    input_read_gwas_columns = [gwas_col_dict['snp'], gwas_col_dict['chrom'],
                               gwas_col_dict['position'],
                               gwas_col_dict['pvalue'], var_id]
    input_read_eqtl_columns = [eqtl_col_dict['snp'], eqtl_col_dict['pvalue'], var_id]

    # creat manhattan plot file
    gwas_preprocessed_pd = create_manhattan_plot_file(gwas_col_dict, gwas_preprocessed_file, var_id, trait_path)

    df_genes_info_path_pd['h_target_rsid'] = ''
    for gene_index, gene_row in df_genes_info_path_pd.iterrows():
        eqtl_file_path = gene_row['eqtl_path']
        gene_id = gene_row['gene_id']
        gene_file_path = os.path.join(trait_path, f'{gene_id}.json')
        if not utils.file_exists(gene_file_path) and (gene_row['chrom'] is not None):
            if eqtl_file_path is None or pd.isnull(eqtl_file_path) or not utils.file_exists(eqtl_file_path):
                logging.info(f'file no exist {eqtl_file_path}')
            else:
                eqtl_file_path_pd = pd.read_table(eqtl_file_path, sep=const.column_spliter, header=0,
                                                  usecols=input_read_eqtl_columns)
                eqtl_file_path_pd.rename(
                    columns={eqtl_col_dict['snp']: 'eqtl_snp',
                             eqtl_col_dict['pvalue']: 'eqtl_pvalue'}, inplace=True)
                merge_pd = pd.merge(left=gwas_preprocessed_pd,
                                    right=eqtl_file_path_pd,
                                    how='inner', on=var_id)
                merge_pd['p'] = merge_pd['eqtl_pvalue'] + merge_pd[gwas_col_dict['pvalue']]
                target_snp_r2 = merge_pd.loc[merge_pd['p'].idxmin()][gwas_col_dict['snp']]
                gene_row['h_target_rsid'] = target_snp_r2
                # mapping r2 info
                r2_pd = sqr.retrieve_ld(chr=int(gene_row['chrom']), rsid=target_snp_r2, population=population,
                                        rsid_col_name=gwas_col_dict['snp'])
                merge_pd = pd.merge(left=merge_pd, right=r2_pd, how='left', on=gwas_col_dict['snp'])

                if merge_pd.shape[0] == 0:
                    logging.info(f'merge_pd have no data {eqtl_file_path}')
                else:
                    merge_pd['r2'].fillna(0, inplace=True)
                    merge_pd['category'] = pd.cut(x=merge_pd['r2'], bins=[0, 0.2, 0.4, 0.6, 0.8, 1, 2],
                                                  labels=['r2_one', 'r2_two', 'r2_thr', 'r2_four', 'r2_five',
                                                          'target'])
                    merge_pd['category'].fillna('r2_one', inplace=True)
                    merge_pd.sort_values(gwas_col_dict['position'], inplace=True)

                    merge_pd = merge_pd.rename(
                        columns={gwas_col_dict['snp']: output_read_columns_mode['snp'],
                                 gwas_col_dict['pvalue']: output_read_columns_mode['pvalue'],
                                 gwas_col_dict['position']: output_read_columns_mode['position'],
                                 gwas_col_dict['chrom']: output_read_columns_mode['chrom'],
                                 })

                    with open(gene_file_path, 'w', encoding='utf-8') as output_file:
                        output_file.write(
                            'callbackGetPlotsReportData(' + json.dumps(merge_pd[
                                [output_read_columns_mode['snp'], 'eqtl_snp', output_read_columns_mode['pvalue'],
                                 'eqtl_pvalue',
                                 output_read_columns_mode['position'], output_read_columns_mode['chrom'], 'r2',
                                 'category',
                                 var_id, ]].to_dict(
                                orient='records'), indent=2) + ')')
                    logging.info(f'gene_file_path : {gene_file_path}')

    with open(os.path.join(trait_path, 'genes_path.json'), 'w', encoding='utf-8') as output_file:
        output_file.write(
            'callbackGetGenesReportData(' + json.dumps(df_genes_info_path_pd.to_dict(
                orient='records'), indent=2) + ')')


def create_manhattan_plot_file(gwas_col_dict=None, gwas_preprocessed_file=None, var_id=None, trait_path=None):
    input_read_gwas_columns = [gwas_col_dict['snp'], gwas_col_dict['chrom'],
                               gwas_col_dict['position'],
                               gwas_col_dict['pvalue'], var_id]

    gwas_preprocessed_pd = \
        pd.read_table(gwas_preprocessed_file, sep=const.column_spliter, header=0,
                      usecols=input_read_gwas_columns)[input_read_gwas_columns]

    df_len = len(gwas_preprocessed_pd)
    manhattan_pf = gwas_preprocessed_pd[
        [gwas_col_dict['snp'], gwas_col_dict['pvalue'], gwas_col_dict['position'], gwas_col_dict['chrom'],
         var_id]]
    manhattan_pf.loc[:, 'category'] = ''
    manhattan_pf.loc[(manhattan_pf[gwas_col_dict['chrom']] % 2 == 0), 'category'] = 'double_chr'
    manhattan_pf.loc[(manhattan_pf[gwas_col_dict['chrom']] % 2 > 0), 'category'] = 'single_chr'
    pval = const.manhattan_min_pval
    pval_max_count = const.manhattan_max_plot
    # utils.clean_data(gwas_df, dup_consider_subset="variant_id")
    if df_len < pval_max_count:
        pval_filter_gwas_df = manhattan_pf.drop(manhattan_pf[manhattan_pf[gwas_col_dict['pvalue']] == 0.0].index,
                                                inplace=False)
    else:
        pval_filter_gwas_df = filter_data_frame_by_p_value(manhattan_pf, pval, pval_max_count, gwas_col_dict['pvalue'],
                                                           inplace=False)
    pval_filter_gwas_df.rename(
        columns={gwas_col_dict['snp']: output_read_columns_mode['snp'],
                 gwas_col_dict['pvalue']: output_read_columns_mode['pvalue'],
                 gwas_col_dict['position']: output_read_columns_mode['position'],
                 gwas_col_dict['chrom']: output_read_columns_mode['chrom']}, inplace=True)
    with open(os.path.join(trait_path, f'data_mant.json'), 'w', encoding='utf-8') as output_file:
        output_file.write(
            'callbackGetManhattanData(' + json.dumps(pval_filter_gwas_df.to_dict(orient='records'), indent=2) + ')')

    return gwas_preprocessed_pd


def filter_data_frame_by_p_value(input_df, p_value, pval_max_count, pvalue_col_name, inplace=True):
    p_value_col_name = pvalue_col_name
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


def __copy_static_files(source_path, target_path):
    if not os.path.exists(target_path):
        os.makedirs(target_path)

    if os.path.exists(source_path):
        # root 所指的是当前正在遍历的这个文件夹的本身的地址
        # dirs 是一个 list，内容是该文件夹中所有的目录的名字(不包括子目录)
        # files 同样是 list, 内容是该文件夹中所有的文件(不包括子目录)
        for root, dirs, files in os.walk(source_path):
            for file in files:
                if len(dirs) > 0:
                    target_file = os.path.join(target_path, root.replace(source_path, ''), file)
                else:
                    target_dir = os.path.join(target_path, root.replace(source_path + '/', ''))
                    if not os.path.exists(target_dir):
                        os.makedirs(target_dir)
                    target_file = os.path.join(target_path, root.replace(source_path + '/', ''), file)

                src_file = os.path.join(root, file)
                shutil.copy(src_file, target_file)


def report_html_handler(output_dir):
    # 复制项目static目录文件到所需展示数据目录下
    current_dir = os.path.dirname(Path(__file__).resolve())
    static_files_path = os.path.join(current_dir, 'static')
    output_index_path = os.path.join(output_dir, 'index.html')
    if not utils.file_exists(output_index_path):
        output_static_files_path = os.path.join(output_dir, 'static')
        __copy_static_files(static_files_path, output_static_files_path)
        # 复制html文件所需展示数据目录下
        shutil.copy(os.path.join(current_dir, 'index.html'), output_index_path)

# if __name__ == '__main__':
# study = const.default_study
# default_config = const.default_config
# config_holder = cf.ConfigHolder(single_config_file=default_config, study=study)
# list_a = [{'trait': 'cell', 'tool_name': 'jlim',
#            'report_path': '/Users/luoyujiao/Downloads/test/jlim_output_20220609110744.tsv',
#            'cfg_pro': config_holder}, {'trait': 'cell2', 'tool_name': 'smr',
#                                    'report_path': '/Users/luoyujiao/Downloads/test/smr_output_20220901182724.tsv',
#                                    'cfg_pro': config_holder}, ]
# report_data_process([])
