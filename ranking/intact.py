import os
from pathlib import Path
import logging
import sys
import pandas as pd

from ranking.constants import GENE_ID_COL_NAME
from ranking.constants import TOOL_SIG_COL_INFO
from ranking.constants import RESULT_TYPE_PVAL
from ranking.constants import AVG_RANKING_COL_NAME


def read_tool_result(rpt, tool_name, sig_col_name):
    
    
    if rpt is None or (not os.path.exists(rpt)) or os.path.getsize(rpt) <= 0:
        return None
    additional_reading_cols = []
    if tool_name == 'predixcan':
        additional_reading_cols = ['zscore']
    if tool_name == 'twas':
        additional_reading_cols = ['TWAS.Z']
    elif tool_name == 'smr':
        additional_reading_cols = ['b_SMR', 'se_SMR', 'p_HEIDI']
    report_df = pd.read_table(rpt, usecols=[GENE_ID_COL_NAME, sig_col_name] + additional_reading_cols)
    report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    if tool_name == 'predixcan':
        report_df.rename(columns={'zscore': f'{tool_name}_zscore'}, inplace=True)
    if tool_name == 'twas':
        report_df.rename(columns={'TWAS.Z': f'{tool_name}_zscore'}, inplace=True)
    elif tool_name == 'smr':
        report_df[f'{tool_name}_zscore'] = report_df['b_SMR'] / report_df['se_SMR']
        report_df.drop(columns=['b_SMR', 'se_SMR'], inplace=True)
    report_df.rename(columns={sig_col_name: tool_name}, inplace=True)
    return report_df
    # result_df has 2 or more columns: [GENE_ID_COL_NAME, tool, tool_zscore], tool col means tool result value


def prepare_ranking_input(output_file_path, rpts):
    
    
    if rpts is None:
        logging.warning('No reports provided')
        return None
    non_empty_rpt_cnt = 0
    for rpt in rpts.values():
        non_empty_rpt_cnt += 0 if rpt is None else 1
    if non_empty_rpt_cnt < 2:
        logging.warning(
            f'At least 2 reports required to run intact, now only got {non_empty_rpt_cnt}, nothing to do')
        return None
    df_list = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        tool_rpt = rpts.get(tool)
        if tool_rpt is None or len(tool_rpt) == 0:
            continue
        df_list.append(read_tool_result(tool_rpt, tool_name=tool, sig_col_name=sig_column))
    ranking_df = None
    for rpt_df in [df for df in df_list if df is not None and df.shape[0] > 0]:
        if ranking_df is None:
            ranking_df = rpt_df
        else:
            ranking_df = pd.merge(left=ranking_df, right=rpt_df,
                                  left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                                  how='outer')
    if ranking_df is None:
        return None
    tool_ranking_cols = []
    prob_tools = []
    pval_tools = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        if tool in ranking_df.columns:
            if sig_type == RESULT_TYPE_PVAL:
                pval_tools.append(tool)
            else:
                prob_tools.append(tool)
            tool_ranking_cols.append(f'{tool}_ranking')
            # 补充各个工具的ranking数据
            ranking_df.sort_values(tool, ascending=sig_type == RESULT_TYPE_PVAL, inplace=True)
            sig_na_bool_series = ranking_df[tool].isna()
            ranking_df.loc[~sig_na_bool_series, f'{tool}_ranking'] = \
                range(1, ranking_df[~sig_na_bool_series].shape[0] + 1)
            ranking_df.loc[sig_na_bool_series, f'{tool}_ranking'] = ranking_df.shape[0] + 1
    ranking_df[AVG_RANKING_COL_NAME] = ranking_df[tool_ranking_cols].mean(axis=1)
    ranking_df.drop(columns=tool_ranking_cols, inplace=True)
    # 如果部分基因在TWAS或者colocalization这两边任何一方没有结果, 怎么处理? input如果有一方是NA, 那么output也是NA
    if len(prob_tools) == 0 or len(pval_tools) == 0:
        logging.warning('To run intact, you need both TWAS results and colocalization results, skipping')
        return None
    zscore_cols = [f'{col}_zscore' for col in pval_tools]
    ranking_df['probability'] = ranking_df[prob_tools].mean(axis=1)
    # ranking_df.loc[ranking_df['probability'].isna(), 'probability'] = 0
    ranking_df['zscore'] = ranking_df[zscore_cols].mean(axis=1)
    # ranking_df.loc[ranking_df['zscore'].isna(), 'zscore'] = 0
    ranking_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


def run_ranking(rpt=None, output_file_path=None, prior_fun=None):
    
    
    ranking_input_file = os.path.join(os.path.dirname(output_file_path), f'pre_{output_file_path.split(os.sep)[-1]}')
    prepare_ranking_input(output_file_path=ranking_input_file, rpts=rpt)
    if not os.path.exists(ranking_input_file) or os.path.getsize(ranking_input_file) <= 0:
        return None
    intact_input_df = pd.read_table(ranking_input_file, nrows=0)
    if 'probability' in intact_input_df.columns and 'zscore' in intact_input_df.columns:
        rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'intact.R')
        os.system(f'Rscript --no-save --no-restore {rscript_path} {ranking_input_file} {output_file_path} {prior_fun}')
    return output_file_path
