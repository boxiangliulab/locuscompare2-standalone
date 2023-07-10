import os
import pandas as pd
import logging
from pathlib import Path
from ranking.constants import GENE_ID_COL_NAME
from ranking.constants import TOOL_SIG_COL_INFO
from ranking.constants import RESULT_TYPE_PVAL

RANK_COL_NAME = 'rank'


def prepare_rra_input(output_file_path, rpts):
    if rpts is None or len(rpts) == 0:
        logging.warning('No reports provided')
        return None
    if not any(rpts.values()):
        logging.warning('All reports are empty or none')
        return None
    df_list = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        tool_rpt = rpts.get(tool)
        if tool_rpt is None or len(tool_rpt) == 0:
            continue
        df_list.append(read_tool_result(tool_rpt, tool_name=tool, sig_col_name=sig_column, sig_type=sig_type))

    ranking_df = None
    for rpt_df in [df for df in df_list if df is not None and df.shape[0] > 0]:
        if ranking_df is None:
            ranking_df = rpt_df
        else:
            ranking_df = pd.merge(left=ranking_df, right=rpt_df,
                                  left_on=RANK_COL_NAME, right_on=RANK_COL_NAME,
                                  how='outer')
    if ranking_df is None:
        return None
    # RRA trait ranking list as rank based,  the order of the genes is used as the ranking
    ranking_df.sort_values(RANK_COL_NAME, ascending=True, inplace=True)
    ranking_df = ranking_df.reindex(
        columns=[RANK_COL_NAME] + [col for col in ranking_df.columns if col not in [RANK_COL_NAME]], copy=False)
    ranking_df.to_csv(output_file_path, sep='\t', header=True, index=False, na_rep='NA')
    return output_file_path


def read_tool_result(rpt, tool_name, sig_col_name, sig_type):
    if rpt is None or (not os.path.exists(rpt)) or os.path.getsize(rpt) <= 0:
        return None
    rpt_df = pd.read_table(rpt, usecols=[sig_col_name, GENE_ID_COL_NAME])
    rpt_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    rpt_df.sort_values(sig_col_name, ascending=sig_type == RESULT_TYPE_PVAL, inplace=True)
    rpt_df.drop(labels=sig_col_name, axis=1, inplace=True)
    rpt_df[tool_name] = rpt_df[GENE_ID_COL_NAME]
    rpt_df.drop(labels=GENE_ID_COL_NAME, axis=1, inplace=True)
    rpt_df[RANK_COL_NAME] = range(1, rpt_df.shape[0] + 1)
    return rpt_df


def run_ranking(output_file_path, rpt=None, sample_size=None, method=None):
    ranking_input_file = os.path.join(os.path.dirname(output_file_path), f'pre_{output_file_path.split(os.sep)[-1]}')
    prepare_rra_input(output_file_path=ranking_input_file, rpts=rpt)
    if not os.path.exists(ranking_input_file) or os.path.getsize(ranking_input_file) <= 0:
        return None
    rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'rra.R')
    os.system(f'Rscript --no-save --no-restore '
              f'{rscript_path} {ranking_input_file} {output_file_path} {sample_size} {method}')
    return output_file_path
