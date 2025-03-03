import pandas as pd
import os
import sys
from pathlib import Path
import shutil
from datetime import datetime
import argparse
import uuid

PHENOTYPE_ID_COL_NAME = 'phenotype_id'



def calc_threshold_for_pval_rpt(rpt, p_val_col_name, work_dir="~/", fdr_thresh=0.05):
    print(f"calc_threshold_for_pval_rpt")
    if rpt is None or len(rpt) == 0:
        print(f"calc_threshold_for_pval_rpt {rpt} is None")
        return 0, "No result found"
    if not os.path.exists(rpt) or os.path.getsize(rpt) <= 0:
        print(f"calc_threshold_for_pval_rpt {rpt} file not exist")
        return 0, "No result found"
    rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'cal_pval_fdr.r')
    # rscript_path = '/home/liulab/codes/locuscompare2-standalone/fdr/cal_pval_fdr.r'
    output_file = os.path.join(str(work_dir), 'qvalue.txt')

    print(f"Using R script at: {rscript_path}")
    print(f"rscript_path: {rscript_path}")
    print(f"rpt: {rpt}")
    print(f"output_file: {output_file}")
    print(f"p_val_col_name: {p_val_col_name}")
    print(f"fdr_thresh: {fdr_thresh}")
    os.system(f'Rscript --no-save --no-restore {rscript_path} {rpt} {output_file} {p_val_col_name} {fdr_thresh}')
    if os.path.exists(output_file):
        thresh_df = pd.read_table(output_file)
        # utils.delete_file_if_exists(output_file)
        fdr = thresh_df.at[0, 'fdr']
        pvalue = thresh_df.at[0, 'pvalue']
        print(f'Threshold for {rpt} at fdr {fdr_thresh}:\nfdr:{fdr}\npvalue:{pvalue}')
        return pvalue, f"Default: {pvalue} (FDR < 0.05)"
    else:
        return 0, "FDR threshold calculation failed"

if __name__ == '__main__':
    start_time = datetime.now()
    print(f'start run all coloctools, start time: {start_time}')
    parser = argparse.ArgumentParser()
    # parser.add_argument('--tools', dest='tools_list', default=['all'], type=str, nargs='*', help='coloc tools')
    # config.yml directory
    parser.add_argument('--rpt', dest='rpt', help="")
    parser.add_argument('--p_val_col_name', dest='p_val_col_name', help="")
    parser.add_argument('--work_dir', dest='work_dir', help="")

    parse_args = parser.parse_args()
    calc_threshold_for_pval_rpt(parse_args.rpt, parse_args.p_val_col_name, 
        parse_args.work_dir)

