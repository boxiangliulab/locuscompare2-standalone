import logging
import os
import uuid
from pathlib import Path
import sys
import pandas as pd

from ranking.constants import GENE_ID_COL_NAME
from ranking.constants import RESULT_TYPE_PVAL
from ranking.constants import TOOL_SIG_COL_INFO
from common import coloc_utils as utils


def calc_threshold_for_pval_rpt(rpt, p_val_col_name, work_dir='', fdr_thresh=0.05):
    
    
    if rpt is None or len(rpt) == 0:
        return 0
    if not os.path.exists(rpt) or os.path.getsize(rpt) <= 0:
        return 0
    rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'calc.R')
    output_file = os.path.join(work_dir, str(uuid.uuid4()))
    os.system(f'Rscript --no-save --no-restore {rscript_path} {rpt} {output_file} {p_val_col_name} {fdr_thresh}')
    thresh_df = pd.read_table(output_file)
    utils.delete_file_if_exists(output_file)
    fdr = thresh_df.at[0, 'fdr']
    pvalue = thresh_df.at[0, 'pvalue']
    logging.info(f'Threshold for {rpt} at fdr {fdr_thresh}:\nfdr:{fdr}\npvalue:{pvalue}')
    return thresh_df.at[0, 'pvalue']


def calc_threshold_for_prob_rpt(rpt, prob_col_name, fdr_thresh=0.05):
    
    
    if rpt is None or len(rpt) == 0:
        return 1
    if not os.path.exists(rpt) or os.path.getsize(rpt) <= 0:
        return 1
    # ecaviar threshold fix to 0.01
    if prob_col_name == 'clpp':
        return 0.01
    df = pd.read_table(rpt, usecols=[GENE_ID_COL_NAME, prob_col_name])
    if df.empty:
        return 1
    df.sort_values(by=[prob_col_name], ascending=False, inplace=True)
    df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    df['cumulative_count'] = range(1, df.shape[0] + 1)
    df['fdr'] = (1 - df[prob_col_name]).cumsum() / df['cumulative_count']
    threshold_idx = (df['fdr'] - fdr_thresh).abs().idxmin()
    fdr = df['fdr'].loc[threshold_idx]
    prob_thresh = df[prob_col_name].loc[threshold_idx]
    logging.info(f'Threshold for {rpt} at fdr {fdr_thresh}:\nfdr:{fdr}\nprobability:{prob_thresh}')
    return prob_thresh


def calc_threshold(rpt_obj=None, work_dir=''):
    
    
    thresholds = {}
    if rpt_obj is None or len(rpt_obj) == 0:
        logging.warning('No reports provided')
        return thresholds
    if not os.path.exists(work_dir):
        Path(work_dir).mkdir(parents=True, exist_ok=True)
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        tool_rpt = rpt_obj.get(tool)
        if tool_rpt is None or len(tool_rpt) == 0:
            logging.warning(f'{tool_rpt} does not exist, skipping {tool}')
            continue
        if not os.path.exists(tool_rpt) or os.path.getsize(tool_rpt) <= 0:
            logging.warning(f'{tool_rpt} does not exist or is empty, skipping {tool}')
            continue
        if sig_type == RESULT_TYPE_PVAL:
            thresholds[tool] = calc_threshold_for_pval_rpt(tool_rpt, sig_column, work_dir=work_dir)
        else:
            thresholds[tool] = calc_threshold_for_prob_rpt(tool_rpt, sig_column)
    return thresholds
