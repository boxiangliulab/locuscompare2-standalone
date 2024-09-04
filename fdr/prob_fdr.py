import pandas as pd
import os
import sys
from pathlib import Path
import shutil
from datetime import datetime
import argparse
import uuid


GENE_ID_COL_NAME = 'gene_id'


def calc_threshold_for_prob_rpt(rpt, prob_col_name, fdr_thresh=0.05):
    print(f"calc_threshold_for_prob_rpt")
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
    print(f'Threshold for {rpt} at fdr {fdr_thresh}:\nfdr:{fdr}\nprobability:{prob_thresh}')
    return prob_thresh

if __name__ == '__main__':
    start_time = datetime.now()
    print(f'start run all coloctools, start time: {start_time}')
    parser = argparse.ArgumentParser()
    parser.add_argument('--rpt', dest='rpt', help="")
    parser.add_argument('--p_val_col_name', dest='p_val_col_name', help="")

    parse_args = parser.parse_args()
    calc_threshold_for_prob_rpt(parse_args.rpt, parse_args.p_val_col_name)

