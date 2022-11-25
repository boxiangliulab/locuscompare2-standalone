import os
import pandas as pd
from pathlib import Path
from ranking.constants import GENE_ID_COL_NAME
from ranking.constants import TOOL_SIG_COL_INFO
from ranking.constants import RESULT_TYPE_PVAL


def prepare_birra_input(output_file_path, rpts):
    if rpts is None or len(rpts) == 0:
        return
    if not any(rpts.values()):
        return
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
                                  left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                                  how='outer')
    # BIRRA trait ranking as p-values, so we set missing values to max p-value(i.e.:1)
    ranking_df.fillna(1, inplace=True)
    ranking_df = ranking_df.reindex(
        columns=[GENE_ID_COL_NAME] + [col for col in ranking_df.columns if col not in [GENE_ID_COL_NAME]], copy=False)
    ranking_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


def read_tool_result(rpt, tool_name, sig_col_name, sig_type):
    if rpt is None:
        return None
    rpt_df = pd.read_table(rpt, usecols=[sig_col_name, GENE_ID_COL_NAME])
    rpt_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    if sig_type == RESULT_TYPE_PVAL:
        rpt_df.rename(columns={sig_col_name: f'{tool_name}'}, inplace=True)
    else:
        # Convert probability to p-value. TODO this method is not good
        rpt_df[tool_name] = 1 - rpt_df[sig_col_name]
        rpt_df.drop(labels=sig_col_name, axis=1, inplace=True)
    return rpt_df


def run_ranking(output_file_path, rpt=None):
    ranking_input_file = os.path.join(os.path.dirname(output_file_path), f'pre_{output_file_path.split(os.sep)[-1]}')
    prepare_birra_input(output_file_path=ranking_input_file, rpts=rpt)
    rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'birra.R')
    os.system(f'Rscript --no-save --no-restore {rscript_path} {ranking_input_file} {output_file_path}')
    return output_file_path
