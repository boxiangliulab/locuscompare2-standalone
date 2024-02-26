import argparse
import os
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve, average_precision_score, precision_recall_curve
from upsetplot import plot, from_indicators
from venn import pseudovenn
import sys
import ranking.intact as intact
import ranking.rra as rra
from ranking.constants import TOOL_SIG_COL_INFO, GENE_ID_COL_NAME, RESULT_TYPE_PVAL, RESULT_TYPE_PROB
from threshold import calc

prob_col_name = 'PROB'
is_positive_col_name = 'IS_POSITIVE'


def retrieve_std_df(generated_list, sec_causal_type):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    generated_list_df = pd.read_table(generated_list,
                                      usecols=['gene', 'gwas_causal_snp', 'gwas_max_assoc_p', 'gwas_r2',
                                               'eqtl_max_assoc_p_1', f'eqtl_max_assoc_p_{sec_causal_type}'])
    positive_df = generated_list_df[(generated_list_df['gwas_max_assoc_p'] <= 1.0E-5) & (
            generated_list_df['gwas_r2'] >= 0.8) & (generated_list_df['eqtl_max_assoc_p_1'] <= 0.01)].copy()
    positive_df[is_positive_col_name] = 1
    positive_df[GENE_ID_COL_NAME] = positive_df['gene'] + '_1'

    negative_df = generated_list_df[(generated_list_df['gwas_max_assoc_p'] <= 1.0E-5) & (
            generated_list_df['gwas_r2'] >= 0.8) & (generated_list_df[
                                                        f'eqtl_max_assoc_p_{sec_causal_type}'] <= 0.01)].copy()
    negative_df[is_positive_col_name] = 0
    negative_df[GENE_ID_COL_NAME] = negative_df['gene'] + f'_{sec_causal_type}'

    std_df = pd.concat(
        [positive_df[[GENE_ID_COL_NAME, is_positive_col_name]], negative_df[[GENE_ID_COL_NAME, is_positive_col_name]]])
    return std_df


def retrieve_positive_df(generated_list):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    generated_list_df = pd.read_table(generated_list,
                                      usecols=['gene', 'gwas_max_assoc_p', 'gwas_r2',
                                               'eqtl_max_assoc_p_1'])
    positive_df = generated_list_df[(generated_list_df['gwas_max_assoc_p'] <= 1.0E-5) & (
            generated_list_df['gwas_r2'] >= 0.8) & (generated_list_df['eqtl_max_assoc_p_1'] <= 0.01)].copy()
    positive_df.rename(columns={'gene': GENE_ID_COL_NAME}, inplace=True)
    positive_df.drop(columns=[col for col in positive_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    positive_df[is_positive_col_name] = 1
    positive_df.reset_index(drop=True, inplace=True)
    return positive_df


def retrieve_negative_df(generated_list, sec_causal_type):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    generated_list_df = pd.read_table(generated_list,
                                      usecols=['gene', 'gwas_causal_snp', 'gwas_max_assoc_p', 'gwas_r2',
                                               f'eqtl_max_assoc_p_{sec_causal_type}'])
    negative_df = generated_list_df[(generated_list_df['gwas_max_assoc_p'] <= 1.0E-5) & (
            generated_list_df['gwas_r2'] >= 0.8) & (generated_list_df[
                                                        f'eqtl_max_assoc_p_{sec_causal_type}'] <= 0.01)].copy()
    negative_df.rename(columns={'gene': GENE_ID_COL_NAME}, inplace=True)
    negative_df.drop(columns=[col for col in negative_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    negative_df[is_positive_col_name] = 0
    negative_df.reset_index(drop=True, inplace=True)
    return negative_df


def prepare_plot_data(generated_list, h1_report, sec_report, sec_causal_type,
                      rpt_prob_col_name=None, rpt_pval_col_name=None, tool=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt_prob_col_name is None and rpt_pval_col_name is None:
        raise ValueError('p-value column name and probability column name can not be both null')
    std_df = retrieve_std_df(generated_list, sec_causal_type)

    if rpt_prob_col_name is None:
        reading_cols = [rpt_pval_col_name, GENE_ID_COL_NAME]
        if tool == 'predixcan':
            reading_cols.append('zscore')
        elif tool == 'twas':
            reading_cols.append('TWAS.Z')
        elif tool == 'smr':
            reading_cols.append('b_SMR')
            reading_cols.append('se_SMR')
            reading_cols.append('p_HEIDI')
        if h1_report is not None and Path(h1_report).exists() and os.path.getsize(h1_report) > 0:
            h1_report_df = pd.read_table(h1_report, usecols=reading_cols)
            h1_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
            h1_report_df[prob_col_name] = 1 - h1_report_df[rpt_pval_col_name]
        else:
            h1_report_df = pd.DataFrame(columns=[rpt_pval_col_name, prob_col_name, GENE_ID_COL_NAME])
        if sec_report is not None and Path(sec_report).exists() and os.path.getsize(sec_report) > 0:
            sec_report_df = pd.read_table(sec_report, usecols=reading_cols)
            sec_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
            # Convert p-value to probability. TODO this method is not good
            sec_report_df[prob_col_name] = 1 - sec_report_df[rpt_pval_col_name]
        else:
            sec_report_df = pd.DataFrame(columns=[rpt_pval_col_name, prob_col_name, GENE_ID_COL_NAME])
    else:
        if h1_report is not None and Path(h1_report).exists() and os.path.getsize(h1_report) > 0:
            h1_report_df = pd.read_table(h1_report, usecols=[rpt_prob_col_name, GENE_ID_COL_NAME])
            h1_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
            h1_report_df[prob_col_name] = h1_report_df[rpt_prob_col_name]
        else:
            h1_report_df = pd.DataFrame(columns=[rpt_prob_col_name, prob_col_name, GENE_ID_COL_NAME])
        if sec_report is not None and Path(sec_report).exists() and os.path.getsize(sec_report) > 0:
            sec_report_df = pd.read_table(sec_report, usecols=[rpt_prob_col_name, GENE_ID_COL_NAME])
            sec_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
            sec_report_df[prob_col_name] = sec_report_df[rpt_prob_col_name]
        else:
            sec_report_df = pd.DataFrame(columns=[rpt_prob_col_name, prob_col_name, GENE_ID_COL_NAME])
    h1_report_df[GENE_ID_COL_NAME] = h1_report_df[GENE_ID_COL_NAME] + '_1'
    sec_report_df[GENE_ID_COL_NAME] = sec_report_df[GENE_ID_COL_NAME] + f'_{sec_causal_type}'
    report_df = pd.concat([h1_report_df, sec_report_df])
    # TODO MHY is it right to use how as outer?
    result_df = pd.merge(left=std_df, right=report_df,
                         left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                         how='outer')
    tp_na_bool_series = result_df[is_positive_col_name].isna()
    result_df.loc[result_df[tp_na_bool_series].index, is_positive_col_name] = 0
    prob_na_bool_series = result_df[prob_col_name].isna()
    result_df.loc[result_df[prob_na_bool_series].index, prob_col_name] = 0
    if rpt_prob_col_name is None:
        # TODO MHY here 这里还在补值, plot value不能为NA
        result_df.loc[result_df[prob_na_bool_series].index, rpt_pval_col_name] = 1
    else:
        result_df.loc[result_df[prob_na_bool_series].index, rpt_prob_col_name] = 0
    # result_df has 4 columns: [gene_id, prob_col_name, is_positive_col_name, (rpt_prob_col_name or rpt_pval_col_name)]
    return result_df


def plot_single_roc(generated_list, h1_report, sec_report, sec_causal_type,
                    rpt_prob_col_name=None, rpt_pval_col_name=None, tool=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    result_df = prepare_plot_data(generated_list, h1_report, sec_report, sec_causal_type, rpt_prob_col_name,
                                  rpt_pval_col_name, tool=tool)
    plot_roc_curve(result_df[is_positive_col_name], result_df[prob_col_name], tool=tool)


def plot_roc_curve(true_y, y_prob, tool):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    auc = roc_auc_score(true_y, y_prob)
    if tool is None:
        lbl = 'AUC = {:0.4f}'.format(auc)
    else:
        lbl = tool + ': AUC = {:0.4f}'
        lbl = lbl.format(auc)
    plt.plot(fpr, tpr, label=lbl)
    plt.legend(loc='lower right')


def plot_all_roc(generated_file_path, h1_rpt_obj, sec_rpt_obj, sec_causal_type=1, output_figure_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    plt.figure().clear()
    xlocs, xlabels = plt.xticks()
    for idx, loc in enumerate(xlocs):
        xlabels[idx] = '{:0.1f}'.format(1 - loc)
    # change xtick labels from FPR to Specificity
    plt.xticks(xlocs, xlabels)
    plt.xlim(-0.05, 1.05)
    plt.ylabel('Sensitivity')
    plt.xlabel('Specificity')
    plt.plot([0, 1], [0, 1], color='grey', linewidth=1, linestyle='-')
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            plot_single_roc(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                            rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                            rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                            tool=tool)
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


def plot_single_prc(generated_list, h1_report, sec_report, sec_causal_type,
                    rpt_prob_col_name=None, rpt_pval_col_name=None, tool=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    result_df = prepare_plot_data(generated_list, h1_report, sec_report, sec_causal_type, rpt_prob_col_name,
                                  rpt_pval_col_name)
    plot_prc_curve(result_df[is_positive_col_name], result_df[prob_col_name], tool=tool)


def plot_all_prc(generated_file_path, h1_rpt_obj, sec_rpt_obj, sec_causal_type=1, output_figure_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    plt.figure().clear()
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            plot_single_prc(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                            rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                            rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                            tool=tool)
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


def plot_prc_curve(true_y, y_prob, tool):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    precision, recall, thresholds = precision_recall_curve(true_y, y_prob)
    ap = average_precision_score(true_y, y_prob)
    if tool is None:
        lbl = 'AP = {:0.4f}'.format(ap)
    else:
        lbl = tool + ': AP = {:0.4f}'
        lbl = lbl.format(ap)
    plt.plot(recall, precision, label=lbl)
    plt.legend(loc='upper right')


def plot_all_against_ensemble_roc(generated_file_path, h1_rpt_obj, sec_rpt_obj, sec_causal_type=1,
                                  output_figure_path=None,
                                  output_dir=''):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    plt.figure().clear()
    xlocs, xlabels = plt.xticks()
    for idx, loc in enumerate(xlocs):
        xlabels[idx] = '{:0.1f}'.format(1 - loc)
    # change xtick labels from FPR to Specificity
    plt.xticks(xlocs, xlabels)
    plt.xlim(-0.05, 1.05)
    plt.ylabel('Sensitivity')
    plt.xlabel('Specificity')
    plt.plot([0, 1], [0, 1], color='grey', linewidth=1, linestyle='-')
    rpt_obj = {}
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            tool_df = prepare_plot_data(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                                        rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                        rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                        tool=tool)
            if tool_df.empty:
                continue
            plot_roc_curve(tool_df[is_positive_col_name], tool_df[prob_col_name], tool=tool)
            tool_output = os.path.join(output_dir, f'{tool}_result.tsv')
            tool_df.to_csv(tool_output, sep='\t', header=True, index=False)
            rpt_obj[tool] = tool_output
    if len(rpt_obj) == 0:
        print('all report is empty, nothing to do')
        return
    # birra.run_ranking(output_file_path=ensemble_output_file, rpt=rpt_obj)
    # ranking_result_df = pd.read_table(ensemble_output_file, usecols=[GENE_ID_COL_NAME, 'birra_ranking'])
    # Convert birra ranking to probability. TODO this method is not good
    # ranking_result_df[prob_col_name] = 1 - ranking_result_df['birra_ranking'] / ranking_result_df['birra_ranking'].max()
    std_df = retrieve_std_df(generated_file_path, sec_causal_type)
    # ranking_result_df = pd.merge(left=ranking_result_df, right=std_df,
    #                              left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                              how='left')
    # tp_na_bool_series = ranking_result_df[is_positive_col_name].isna()
    # ranking_result_df.loc[ranking_result_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(ranking_result_df[is_positive_col_name], ranking_result_df[prob_col_name], tool='birra_ranking')
    # --------
    rra_geo_output_file = os.path.join(output_dir, 'rra_geo.tsv')
    rra.run_ranking(output_file_path=rra_geo_output_file, rpt=rpt_obj, method='GEO')
    rgeo_df = pd.read_table(rra_geo_output_file, usecols=[GENE_ID_COL_NAME, 'geo_ranking'])
    # Convert rra p_value to probability. TODO this method is not good
    rgeo_df[prob_col_name] = 1 - rgeo_df['geo_ranking']
    rgeo_df = pd.merge(left=rgeo_df, right=std_df,
                       left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                       how='left')
    tp_na_bool_series = rgeo_df[is_positive_col_name].isna()
    rgeo_df.loc[rgeo_df[tp_na_bool_series].index, is_positive_col_name] = 0
    plot_roc_curve(rgeo_df[is_positive_col_name], rgeo_df[prob_col_name], tool='rGEO')
    # --------
    # rra_stuart_output_file = os.path.join(output_dir, 'rra_stuart.tsv')
    # rra.run_ranking(output_file_path=rra_stuart_output_file, rpt=rpt_obj, sample_size=191647, method='stuart')
    # stuart_df = pd.read_table(rra_stuart_output_file, usecols=[GENE_ID_COL_NAME, 'stuart_p_value'])
    # # Convert rra p_value to probability. TODO this method is not good
    # stuart_df[prob_col_name] = 1 - stuart_df['stuart_p_value']
    # stuart_df = pd.merge(left=stuart_df, right=std_df,
    #                      left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                      how='left')
    # tp_na_bool_series = stuart_df[is_positive_col_name].isna()
    # stuart_df.loc[stuart_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(stuart_df[is_positive_col_name], stuart_df[prob_col_name], tool='stuart')
    # --------xgboost lr
    # train_output_file = os.path.join(output_dir, 'train_xgb_lr.tsv')
    # machine_learning.run_ranking(rpt=rpt_obj, output_file_path=train_output_file)
    # train_df = pd.read_table(train_output_file, usecols=[GENE_ID_COL_NAME, 'ml_probability'])
    # train_df = pd.merge(left=train_df, right=std_df,
    #                     left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                     how='left')
    # tp_na_bool_series = train_df[is_positive_col_name].isna()
    # train_df.loc[train_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_df[is_positive_col_name], train_df['ml_probability'], tool='ensemble')
    # --------xgboost rf
    # train_xgb_rf_output_file = os.path.join(output_dir, 'train_xgb_rf.tsv')
    # machine_learning.run_ranking_xgb_rf(rpt=rpt_obj, output_file_path=train_xgb_rf_output_file)
    # train_xgb_rf_df = pd.read_table(train_xgb_rf_output_file, usecols=[GENE_ID_COL_NAME, 'ml_probability'])
    # train_xgb_rf_df = pd.merge(left=train_xgb_rf_df, right=std_df,
    #                            left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                            how='left')
    # tp_na_bool_series = train_xgb_rf_df[is_positive_col_name].isna()
    # train_xgb_rf_df.loc[train_xgb_rf_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_xgb_rf_df[is_positive_col_name], train_xgb_rf_df['ml_probability'], tool='ensemble xgb rf')
    # --------scikit-learn svm
    # train_svm_output_file = os.path.join(output_dir, 'train_svm.tsv')
    # machine_learning.run_ranking_svm(rpt=rpt_obj, output_file_path=train_svm_output_file)
    # train_svm_df = pd.read_table(train_svm_output_file, usecols=[GENE_ID_COL_NAME, 'ml_probability'])
    # train_svm_df = pd.merge(left=train_svm_df, right=std_df,
    #                         left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                         how='left')
    # tp_na_bool_series = train_svm_df[is_positive_col_name].isna()
    # train_svm_df.loc[train_svm_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_svm_df[is_positive_col_name], train_svm_df['ml_probability'], tool='ensemble svm')
    # --------scikit-learn random forest
    # train_rf_output_file = os.path.join(output_dir, 'train_rf.tsv')
    # machine_learning.run_ranking_rf(rpt=rpt_obj, output_file_path=train_rf_output_file)
    # train_rf_df = pd.read_table(train_rf_output_file, usecols=[GENE_ID_COL_NAME, 'ml_probability'])
    # train_rf_df = pd.merge(left=train_rf_df, right=std_df,
    #                        left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                        how='left')
    # tp_na_bool_series = train_rf_df[is_positive_col_name].isna()
    # train_rf_df.loc[train_rf_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_rf_df[is_positive_col_name], train_rf_df['ml_probability'], tool='ensemble rf')
    # --------scikit-learn logistic regression
    # train_lr_output_file = os.path.join(output_dir, 'train_lr.tsv')
    # machine_learning.run_ranking_lr(rpt=rpt_obj, output_file_path=train_lr_output_file)
    # train_lr_df = pd.read_table(train_lr_output_file, usecols=[GENE_ID_COL_NAME, 'ml_probability'])
    # train_lr_df = pd.merge(left=train_lr_df, right=std_df,
    #                        left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                        how='left')
    # tp_na_bool_series = train_lr_df[is_positive_col_name].isna()
    # train_lr_df.loc[train_lr_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_lr_df[is_positive_col_name], train_lr_df['ml_probability'], tool='ensemble lr')
    # --------intact
    # intact_output_file = os.path.join(output_dir, 'intact_linear.tsv')
    # intact.run_ranking(rpt=rpt_obj, output_file_path=intact_output_file, prior_fun='linear')
    # intact_df = pd.read_table(intact_output_file, usecols=[GENE_ID_COL_NAME, 'intact_probability'])
    # intact_df = pd.merge(left=intact_df, right=std_df,
    #                      left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                      how='left')
    # tp_na_bool_series = intact_df[is_positive_col_name].isna()
    # intact_df.loc[intact_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # intact_df.loc[intact_df['intact_probability'].isna(), 'intact_probability'] = 0
    # plot_roc_curve(intact_df[is_positive_col_name], intact_df['intact_probability'], tool='intact linear')
    # --------intact_hybrid
    # intact_hybrid_output_file = os.path.join(output_dir, 'intact_hybrid.tsv')
    # intact.run_ranking(rpt=rpt_obj, output_file_path=intact_hybrid_output_file, prior_fun='hybrid')
    # intact_hybrid_df = pd.read_table(intact_hybrid_output_file, usecols=[GENE_ID_COL_NAME, 'intact_probability'])
    # intact_hybrid_df = pd.merge(left=intact_hybrid_df, right=std_df,
    #                             left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                             how='left')
    # tp_na_bool_series = intact_hybrid_df[is_positive_col_name].isna()
    # intact_hybrid_df.loc[intact_hybrid_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # intact_hybrid_df.loc[intact_hybrid_df['intact_probability'].isna(), 'intact_probability'] = 0
    # plot_roc_curve(intact_hybrid_df[is_positive_col_name], intact_hybrid_df['intact_probability'], tool='intact')
    # --------intact_step
    # intact_step_output_file = os.path.join(output_dir, 'intact_step.tsv')
    # intact.run_ranking(rpt=rpt_obj, output_file_path=intact_step_output_file, prior_fun='step')
    # intact_step_df = pd.read_table(intact_step_output_file, usecols=[GENE_ID_COL_NAME, 'intact_probability'])
    # intact_step_df = pd.merge(left=intact_step_df, right=std_df,
    #                           left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                           how='left')
    # tp_na_bool_series = intact_step_df[is_positive_col_name].isna()
    # intact_step_df.loc[intact_step_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # intact_step_df.loc[intact_step_df['intact_probability'].isna(), 'intact_probability'] = 0
    # plot_roc_curve(intact_step_df[is_positive_col_name], intact_step_df['intact_probability'], tool='intact step')
    # --------intact_expit
    intact_expit_output_file = os.path.join(output_dir, 'intact_expit.tsv')
    intact.run_ranking(rpt=rpt_obj, output_file_path=intact_expit_output_file, prior_fun='expit')
    intact_expit_df = pd.read_table(intact_expit_output_file, usecols=[GENE_ID_COL_NAME, 'intact_probability'])
    intact_expit_df = pd.merge(left=intact_expit_df, right=std_df,
                               left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                               how='left')
    tp_na_bool_series = intact_expit_df[is_positive_col_name].isna()
    intact_expit_df.loc[intact_expit_df[tp_na_bool_series].index, is_positive_col_name] = 0
    intact_expit_df.loc[intact_expit_df['intact_probability'].isna(), 'intact_probability'] = 0
    plot_roc_curve(intact_expit_df[is_positive_col_name], intact_expit_df['intact_probability'], tool='intact expit')
    # ---------simple ranking
    # simple_ranking_result_df = pd.read_table(intact_expit_output_file, usecols=[GENE_ID_COL_NAME, 'avg_ranking'])
    # # Convert simple mean ranking to probability. TODO this method is not good
    # simple_ranking_result_df[prob_col_name] = 1 - simple_ranking_result_df['avg_ranking'] / simple_ranking_result_df[
    #     'avg_ranking'].max()
    # std_df = retrieve_std_df(generated_file_path, sec_causal_type)
    # simple_ranking_result_df = pd.merge(left=simple_ranking_result_df, right=std_df,
    #                                     left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                                     how='left')
    # tp_na_bool_series = simple_ranking_result_df[is_positive_col_name].isna()
    # simple_ranking_result_df.loc[simple_ranking_result_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(simple_ranking_result_df[is_positive_col_name], simple_ranking_result_df[prob_col_name],
    #                tool='avg ranking')
    # --------
    # --------pu learning
    # pu_learning_output_file = os.path.join(output_dir, 'train_pu.tsv')
    # machine_learning.run_ranking_pu(rpt=rpt_obj, output_file_path=pu_learning_output_file)
    # pu_learning_df = pd.read_table(pu_learning_output_file, usecols=[GENE_ID_COL_NAME, 'ml_probability'])
    # pu_learning_df = pd.merge(left=pu_learning_df, right=std_df,
    #                           left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                           how='left')
    # tp_na_bool_series = pu_learning_df[is_positive_col_name].isna()
    # pu_learning_df.loc[pu_learning_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(pu_learning_df[is_positive_col_name], pu_learning_df['ml_probability'], tool='pu learning')
    # --------
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


def plot_all_against_ensemble_prc(generated_file_path, h1_rpt_obj, sec_rpt_obj, sec_causal_type=1,
                                  output_figure_path=None,
                                  output_dir=''):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    plt.figure().clear()
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    rpt_obj = {}
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            tool_df = prepare_plot_data(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                                        rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                        rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                        tool=tool)
            if tool_df.empty:
                continue
            plot_prc_curve(tool_df[is_positive_col_name], tool_df[prob_col_name], tool=tool)
            tool_output = os.path.join(output_dir, f'{tool}_result.tsv')
            tool_df.to_csv(tool_output, sep='\t', header=True, index=False)
            rpt_obj[tool] = tool_output
    if len(rpt_obj) == 0:
        print('all report is empty, nothing to do')
        return
    std_df = retrieve_std_df(generated_file_path, sec_causal_type)
    # --------
    # rra_geo_output_file = os.path.join(output_dir, 'rra_geo.tsv')
    # rra.run_ranking(output_file_path=rra_geo_output_file, rpt=rpt_obj, method='GEO')
    # rgeo_df = pd.read_table(rra_geo_output_file, usecols=[GENE_ID_COL_NAME, 'geo_ranking'])
    # # Convert rra p_value to probability. TODO this method is not good
    # rgeo_df[prob_col_name] = 1 - rgeo_df['geo_ranking']
    # rgeo_df = pd.merge(left=rgeo_df, right=std_df,
    #                    left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                    how='left')
    # tp_na_bool_series = rgeo_df[is_positive_col_name].isna()
    # rgeo_df.loc[rgeo_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_prc_curve(rgeo_df[is_positive_col_name], rgeo_df[prob_col_name], tool='rGEO')
    # --------intact_expit
    intact_expit_output_file = os.path.join(output_dir, 'intact_expit.tsv')
    intact.run_ranking(rpt=rpt_obj, output_file_path=intact_expit_output_file, prior_fun='expit')
    intact_expit_df = pd.read_table(intact_expit_output_file, usecols=[GENE_ID_COL_NAME, 'intact_probability'])
    intact_expit_df = pd.merge(left=intact_expit_df, right=std_df,
                               left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                               how='left')
    tp_na_bool_series = intact_expit_df[is_positive_col_name].isna()
    intact_expit_df.loc[intact_expit_df[tp_na_bool_series].index, is_positive_col_name] = 0
    intact_expit_df.loc[intact_expit_df['intact_probability'].isna(), 'intact_probability'] = 0
    plot_prc_curve(intact_expit_df[is_positive_col_name], intact_expit_df['intact_probability'], tool='intact expit')
    # --------
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


def plot_bar(generated_file_path, h1_rpt_obj, sec_rpt_obj, sec_causal_type=1, output_figure_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    df_list = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            coloc_df = prepare_plot_data(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                                         rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                         rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                         tool=tool)
            if coloc_df.empty:
                continue
            df_list.append(coloc_df)
    if len(df_list) == 0:
        print('all report is empty, nothing to do')
        return
    tool_thresholds = calc.calc_threshold(rpt_obj=h1_rpt_obj, work_dir=os.path.dirname(output_figure_path))
    ranking_df = None
    for rpt_df in df_list:
        rpt_df.drop(columns=[is_positive_col_name, prob_col_name], inplace=True)
        if ranking_df is None:
            ranking_df = rpt_df
        else:
            ranking_df = pd.merge(left=ranking_df, right=rpt_df,
                                  left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                                  how='outer')
    std_df = retrieve_std_df(generated_file_path, sec_causal_type)
    ranking_df = pd.merge(left=std_df, right=ranking_df,
                          left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                          how='outer')
    tool_positive_cols = []
    tool_sensitivity = []
    tool_specificity = []
    tool_f1 = []
    tools = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        if sig_column not in ranking_df.columns:
            continue
        ranking_df.rename(columns={sig_column: tool}, inplace=True)
        # thresholds are covert to tool specific!
        if tool == 'smr':
            ranking_df[f'{tool}_positive'] = (ranking_df[tool] < tool_thresholds[tool]) & (
                    ranking_df['p_HEIDI'] > 0.05)
        elif sig_type == RESULT_TYPE_PROB:
            ranking_df[f'{tool}_positive'] = ranking_df[tool] > tool_thresholds[tool]
        elif sig_type == RESULT_TYPE_PVAL:
            ranking_df[f'{tool}_positive'] = ranking_df[tool] < tool_thresholds[tool]
        else:
            raise ValueError(f'Tool {tool} is not recognized, what is its threshold for it?')
        tp = ranking_df[ranking_df[f'{tool}_positive'] & (ranking_df[is_positive_col_name] == 1)].shape[0]
        # fn = ranking_df[~ranking_df[f'{tool}_positive'] & (ranking_df[is_positive_col_name] == 1)].shape[0]
        p = ranking_df[ranking_df[is_positive_col_name] == 1].shape[0]
        # TPR = tp/p = tp/(tp+fn)
        tn = ranking_df[(~ranking_df[f'{tool}_positive']) & (ranking_df[is_positive_col_name] == 0)].shape[0]
        n = ranking_df[ranking_df[is_positive_col_name] == 0].shape[0]
        fp = ranking_df[(ranking_df[is_positive_col_name] == 0) & ranking_df[f'{tool}_positive']].shape[0]
        fn = ranking_df[(ranking_df[is_positive_col_name] == 1) & (~ranking_df[f'{tool}_positive'])].shape[0]
        tool_sensitivity.append(tp / p)
        # TNR = tn/n = tn/(tn+fp)
        tool_specificity.append(tn / n)
        positive_series = ranking_df[f'{tool}_positive']
        negative_series = ~ positive_series
        ranking_df[f'{tool}_positive'].mask(positive_series, 1, inplace=True)
        ranking_df[f'{tool}_positive'].mask(negative_series, 0, inplace=True)
        tool_positive_cols.append(f'{tool}_positive')
        tool_f1.append(2 * tp / (2 * tp + fp + fn))
        tools.append(f'{tool}\nT={"{:.1e}".format(tool_thresholds[tool])}')
    ranking_df['mvote_positive'] = ranking_df[tool_positive_cols].sum(axis=1) >= 3
    bar_plot_result = os.path.join(os.path.dirname(output_figure_path), f'bar_result.tsv')
    ranking_df.to_csv(bar_plot_result, sep='\t', header=True, index=False, na_rep='NA')
    mv_tp = ranking_df[ranking_df['mvote_positive'] & (ranking_df[is_positive_col_name] == 1)].shape[0]
    mv_p = ranking_df[ranking_df[is_positive_col_name] == 1].shape[0]
    mv_tn = ranking_df[(~ ranking_df['mvote_positive']) & (ranking_df[is_positive_col_name] == 0)].shape[0]
    mv_n = ranking_df[ranking_df[is_positive_col_name] == 0].shape[0]
    mv_fp = ranking_df[(ranking_df[is_positive_col_name] == 0) & ranking_df['mvote_positive']].shape[0]
    mv_fn = ranking_df[(ranking_df[is_positive_col_name] == 1) & (~ranking_df['mvote_positive'])].shape[0]
    mv_f1 = 2 * mv_tp / (2 * mv_tp + mv_fp + mv_fn)
    tools.append(f'm_vote\nT=cnt>=3')
    # tools.append('majority_vote')
    tool_sensitivity.append(mv_tp / mv_p)
    tool_specificity.append(mv_tn / mv_n)
    tool_f1.append(mv_f1)
    print(f'Thresholds:  {tool_thresholds}')
    print(f'Tools: {tools}\nSensitivities: {tool_sensitivity}')
    print(f'Specificities: {tool_specificity}')
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(tools))
    bar_width = 0.25
    # rects = ax.bar(x - bar_width / 2, tool_sensitivity, bar_width, label='Sensitivity')
    # ax.bar_label(rects, fmt='{:0.3f}')
    # rects = ax.bar(x + bar_width / 2, tool_specificity, bar_width, label='Specificity')
    # ax.bar_label(rects, fmt='{:0.3f}')
    rects = ax.bar((x - bar_width), tool_sensitivity, bar_width, label='Sensitivity')
    ax.bar_label(rects, fmt='{:0.2f}')
    rects = ax.bar(x, tool_specificity, bar_width, label='Specificity')
    ax.bar_label(rects, fmt='{:0.2f}')
    rects = ax.bar(x + bar_width, tool_f1, bar_width, label='F1')
    ax.bar_label(rects, fmt='{:0.2f}')
    ax.set_ylabel('Sensitivity/Specificity/F1')
    ax.set_title('Sensitivity/Specificity/F1 of different tools')
    ax.set_xticks(x, tools)
    ax.legend(loc='upper left', ncols=3)
    ax.set_ylim(0, 1.2)

    # -----rendering thresholds start
    # tool_count = len(tools)
    # ax.set_xlim(-bar_width - 0.2, tool_count + 4)
    # for idx, tool in enumerate(tools):
    #     ax.text(tool_count - bar_width, (tool_count - idx) * 0.1,
    #             f'{tool}:{tool_thresholds[tool] if tool in tool_thresholds.keys() else "count>=3"}')
    # -----rendering thresholds end

    # -----plot without label
    # bar_container = ax.bar(tools, tool_sensitivity)
    # ax.set(ylabel='Sensitivity', xlabel='Colocalization tools', title='Sensitivity of different tools', ylim=(0, 1.1))
    # plt.bar(tools, tool_sensitivity)
    # plt.xlabel("Colocalization tools")
    # plt.ylabel("Sensitivity")
    # plt.title("Sensitivity/Specificity of different tools")
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    plt.close()


def calc_max_tpr_minus_fpr_threshold(generated_file_path=None, h1_rpt_obj=None, sec_rpt_obj=None, sec_causal_type=1):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    tool_thresholds = {}
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        tool_df = prepare_plot_data(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                                    rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                    rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                    tool=tool)
        fpr, tpr, thresholds = roc_curve(tool_df[is_positive_col_name], tool_df[prob_col_name])
        prob_thresh = pd.Series(thresholds).loc[pd.Series(tpr - fpr).idxmax()]
        tool_thresholds[tool] = prob_thresh if sig_type == RESULT_TYPE_PROB else 1 - prob_thresh
    return tool_thresholds


def calc_prc_pct_threshold(generated_file_path=None,
                           h1_rpt_obj=None, sec_rpt_obj=None, sec_causal_type=1, prc_pct=0.95):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    tool_thresholds = {}
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        tool_df = prepare_plot_data(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                                    rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                    rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                    tool=tool)
        precision, recall, thresholds = precision_recall_curve(tool_df[is_positive_col_name], tool_df[prob_col_name])
        precision = precision[:-1].copy()
        prob_thresh = pd.Series(thresholds).loc[pd.Series(precision - prc_pct).abs().idxmin()]
        tool_thresholds[tool] = prob_thresh if sig_type == RESULT_TYPE_PROB else 1 - prob_thresh
    return tool_thresholds


def plot_single_venn(
        generated_file_path,
        tool_thresholds=None,
        coloc_rpt=None,
        smr_rpt=None,
        fastenloc_rpt=None,
        predixcan_rpt=None,
        ecaviar_rpt=None,
        twas_rpt=None,
        out_data_prefix=None,
        output_figure_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if coloc_rpt is not None and Path(coloc_rpt).exists() and os.path.getsize(coloc_rpt) > 0:
        coloc_df = pd.read_table(coloc_rpt, usecols=[GENE_ID_COL_NAME, 'overall_H4'])
        coloc_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        coloc_df.drop(labels=coloc_df[coloc_df['overall_H4'] <= tool_thresholds['coloc']].index, inplace=True)
        coloc_df.drop(columns=[col for col in coloc_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        coloc_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])
    if smr_rpt is not None and Path(smr_rpt).exists() and os.path.getsize(smr_rpt) > 0:
        smr_df = pd.read_table(smr_rpt, usecols=[GENE_ID_COL_NAME, 'p_SMR', 'p_HEIDI'])
        smr_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        smr_df.drop(labels=smr_df[(smr_df['p_SMR'] >= tool_thresholds['smr']) | (smr_df['p_HEIDI'] <= 0.05)].index,
                    inplace=True)
        smr_df.drop(columns=[col for col in smr_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        smr_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])
    if fastenloc_rpt is not None and Path(fastenloc_rpt).exists() and os.path.getsize(fastenloc_rpt) > 0:
        fastenloc_df = pd.read_table(fastenloc_rpt, usecols=[GENE_ID_COL_NAME, 'LCP'])
        fastenloc_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        fastenloc_df.drop(labels=fastenloc_df[fastenloc_df['LCP'] <= tool_thresholds['fastenloc']].index, inplace=True)
        fastenloc_df.drop(columns=[col for col in fastenloc_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        fastenloc_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])

    if predixcan_rpt is not None and Path(predixcan_rpt).exists() and os.path.getsize(predixcan_rpt) > 0:
        predixcan_df = pd.read_table(predixcan_rpt, usecols=[GENE_ID_COL_NAME, 'pvalue'])
        predixcan_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        predixcan_df.drop(labels=predixcan_df[predixcan_df['pvalue'] >= tool_thresholds['predixcan']].index,
                          inplace=True)
        predixcan_df.drop(columns=[col for col in predixcan_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        predixcan_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])

    if ecaviar_rpt is not None and Path(ecaviar_rpt).exists() and os.path.getsize(ecaviar_rpt) > 0:
        ecaviar_df = pd.read_table(ecaviar_rpt, usecols=[GENE_ID_COL_NAME, 'clpp'])
        ecaviar_df.sort_values(by='clpp', ascending=False, inplace=True)
        ecaviar_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        ecaviar_df.drop(labels=ecaviar_df[ecaviar_df['clpp'] <= tool_thresholds['ecaviar']].index, inplace=True)
        ecaviar_df.drop(columns=[col for col in ecaviar_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        ecaviar_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])

    if twas_rpt is not None and Path(twas_rpt).exists() and os.path.getsize(twas_rpt) > 0:
        twas_df = pd.read_table(twas_rpt, usecols=[GENE_ID_COL_NAME, 'TWAS.P'])
        twas_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        twas_df.drop(labels=twas_df[twas_df['TWAS.P'] >= tool_thresholds['twas']].index, inplace=True)
        twas_df.drop(columns=[col for col in twas_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        twas_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])

    dataset_dict = {
        f'coloc:{"{:.1e}".format(tool_thresholds["coloc"])}': set(coloc_df[GENE_ID_COL_NAME].values),
        f'smr:{"{:.1e}".format(tool_thresholds["smr"])}': set(smr_df[GENE_ID_COL_NAME].values),
        f'fastenloc:{"{:.1e}".format(tool_thresholds["fastenloc"])}': set(fastenloc_df[GENE_ID_COL_NAME].values),
        f'predixcan:{"{:.1e}".format(tool_thresholds["predixcan"])}': set(predixcan_df[GENE_ID_COL_NAME].values),
        f'ecaviar:{"{:.1e}".format(tool_thresholds["ecaviar"])}': set(ecaviar_df[GENE_ID_COL_NAME].values),
        f'twas:{"{:.1e}".format(tool_thresholds["twas"])}': set(twas_df[GENE_ID_COL_NAME].values)
    }
    pseudovenn(dataset_dict, figsize=(12, 12), hint_hidden=False, cmap="plasma")
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    plt.close()
    # save data to file
    result_df = __merge_tool_gene_id(generated_file_path, coloc_df, smr_df, fastenloc_df, predixcan_df, ecaviar_df,
                                     twas_df)
    pre_venn = os.path.join(os.path.dirname(output_figure_path), f'{out_data_prefix}_pre_venn.tsv')
    result_df.to_csv(pre_venn, sep='\t', header=True, index=False, na_rep='NA')


def plot_venn(generated_file_path=None,
              h1_rpt_obj=None, sec_rpt_obj=None,
              sec_causal_type=1,
              genetic_model=None,
              h1_figure_path=None,
              hgm_figure_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    h1_tool_thresholds = calc.calc_threshold(rpt_obj=h1_rpt_obj, work_dir=os.path.dirname(h1_figure_path))
    sec_tool_thresholds = calc.calc_threshold(rpt_obj=sec_rpt_obj, work_dir=os.path.dirname(h1_figure_path))
    print(f'H1 Thresholds: {h1_tool_thresholds}')
    print(f'SEC Thresholds: {sec_tool_thresholds}')
    plot_single_venn(generated_file_path,
                     h1_tool_thresholds,
                     coloc_rpt=h1_rpt_obj.get('coloc'),
                     smr_rpt=h1_rpt_obj.get('smr'),
                     fastenloc_rpt=h1_rpt_obj.get('fastenloc'),
                     predixcan_rpt=h1_rpt_obj.get('predixcan'),
                     ecaviar_rpt=h1_rpt_obj.get('ecaviar'),
                     twas_rpt=h1_rpt_obj.get('twas'),
                     out_data_prefix='h1',
                     output_figure_path=h1_figure_path)
    plot_single_venn(generated_file_path,
                     sec_tool_thresholds,
                     coloc_rpt=sec_rpt_obj.get('coloc'),
                     smr_rpt=sec_rpt_obj.get('smr'),
                     fastenloc_rpt=sec_rpt_obj.get('fastenloc'),
                     predixcan_rpt=sec_rpt_obj.get('predixcan'),
                     ecaviar_rpt=sec_rpt_obj.get('ecaviar'),
                     twas_rpt=sec_rpt_obj.get('twas'),
                     out_data_prefix=genetic_model,
                     output_figure_path=hgm_figure_path)


def __merge_tool_gene_id(generated_file_path, coloc_df, smr_df, fastenloc_df, predixcan_df, ecaviar_df, twas_df):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    merge_key = '__key'
    coloc_df[merge_key] = coloc_df[GENE_ID_COL_NAME]
    smr_df[merge_key] = smr_df[GENE_ID_COL_NAME]
    fastenloc_df[merge_key] = fastenloc_df[GENE_ID_COL_NAME]
    predixcan_df[merge_key] = predixcan_df[GENE_ID_COL_NAME]
    ecaviar_df[merge_key] = ecaviar_df[GENE_ID_COL_NAME]
    twas_df[merge_key] = twas_df[GENE_ID_COL_NAME]

    positive_df = retrieve_positive_df(generated_file_path)
    positive_df.drop(columns=is_positive_col_name, inplace=True)
    positive_df[merge_key] = positive_df[GENE_ID_COL_NAME]

    result_df = pd.merge(left=positive_df, right=coloc_df,
                         left_on=merge_key, right_on=merge_key,
                         suffixes=(None, "_coloc"),
                         how='left')
    result_df = pd.merge(left=result_df, right=smr_df,
                         left_on=merge_key, right_on=merge_key,
                         suffixes=(None, "_smr"),
                         how='left')
    result_df = pd.merge(left=result_df, right=fastenloc_df,
                         left_on=merge_key, right_on=merge_key,
                         suffixes=(None, "_fastenloc"),
                         how='left')
    result_df = pd.merge(left=result_df, right=predixcan_df,
                         left_on=merge_key, right_on=merge_key,
                         suffixes=(None, "_predixcan"),
                         how='left')
    result_df = pd.merge(left=result_df, right=ecaviar_df,
                         left_on=merge_key, right_on=merge_key,
                         suffixes=(None, "_ecaviar"),
                         how='left')
    result_df = pd.merge(left=result_df, right=twas_df,
                         left_on=merge_key, right_on=merge_key,
                         suffixes=(None, "_twas"),
                         how='left')
    result_df.rename(columns={
        'gene_id_coloc': 'coloc',
        'gene_id_smr': 'smr',
        'gene_id_fastenloc': 'fastenloc',
        'gene_id_predixcan': 'predixcan',
        'gene_id_ecaviar': 'ecaviar',
        'gene_id_twas': 'twas'},
        inplace=True)
    result_df.drop(columns=merge_key, inplace=True)
    return result_df


def plot_mean_sd_combinations_bar(
        generated_file_path=None,
        h1_rpt_obj=None, sec_rpt_obj=None,
        sec_causal_type=1,
        output_figure_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    tool_thresholds = calc.calc_threshold(rpt_obj=h1_rpt_obj, work_dir=os.path.dirname(output_figure_path))
    print(f'Thresholds: {tool_thresholds}')
    h1_coloc_rpt = h1_rpt_obj.get('coloc')
    if h1_coloc_rpt is not None and Path(h1_coloc_rpt).exists() and os.path.getsize(h1_coloc_rpt) > 0:
        coloc_df = pd.read_table(h1_coloc_rpt, usecols=[GENE_ID_COL_NAME, 'overall_H4'])
        coloc_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        coloc_df.drop(labels=coloc_df[coloc_df['overall_H4'] <= tool_thresholds['coloc']].index, inplace=True)
        coloc_df.drop(columns=[col for col in coloc_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        coloc_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])
    h1_smr_rpt = h1_rpt_obj.get('smr')
    if h1_smr_rpt is not None and Path(h1_smr_rpt).exists() and os.path.getsize(h1_smr_rpt) > 0:
        smr_df = pd.read_table(h1_smr_rpt, usecols=[GENE_ID_COL_NAME, 'p_SMR', 'p_HEIDI'])
        smr_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        smr_df.drop(labels=smr_df[(smr_df['p_SMR'] >= tool_thresholds['smr']) | (smr_df['p_HEIDI'] <= 0.05)].index,
                    inplace=True)
        smr_df.drop(columns=[col for col in smr_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        smr_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])
    h1_fastenloc_rpt = h1_rpt_obj.get('fastenloc')
    if h1_fastenloc_rpt is not None and Path(h1_fastenloc_rpt).exists() and os.path.getsize(h1_fastenloc_rpt) > 0:
        fastenloc_df = pd.read_table(h1_fastenloc_rpt, usecols=[GENE_ID_COL_NAME, 'LCP'])
        fastenloc_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        fastenloc_df.drop(labels=fastenloc_df[fastenloc_df['LCP'] <= tool_thresholds['fastenloc']].index, inplace=True)
        fastenloc_df.drop(columns=[col for col in fastenloc_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        fastenloc_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])
    h1_predixcan_rpt = h1_rpt_obj.get('predixcan')
    if h1_predixcan_rpt is not None and Path(h1_predixcan_rpt).exists() and os.path.getsize(h1_predixcan_rpt) > 0:
        predixcan_df = pd.read_table(h1_predixcan_rpt, usecols=[GENE_ID_COL_NAME, 'pvalue'])
        predixcan_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        predixcan_df.drop(labels=predixcan_df[predixcan_df['pvalue'] >= tool_thresholds['predixcan']].index,
                          inplace=True)
        predixcan_df.drop(columns=[col for col in predixcan_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        predixcan_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])
    h1_ecaviar_rpt = h1_rpt_obj.get('ecaviar')
    if h1_ecaviar_rpt is not None and Path(h1_ecaviar_rpt).exists() and os.path.getsize(h1_ecaviar_rpt) > 0:
        ecaviar_df = pd.read_table(h1_ecaviar_rpt, usecols=[GENE_ID_COL_NAME, 'clpp'])
        ecaviar_df.sort_values(by='clpp', ascending=False, inplace=True)
        ecaviar_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        ecaviar_df.drop(labels=ecaviar_df[ecaviar_df['clpp'] <= tool_thresholds['ecaviar']].index, inplace=True)
        ecaviar_df.drop(columns=[col for col in ecaviar_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        ecaviar_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])
    h1_twas_rpt = h1_rpt_obj.get('twas')
    if h1_twas_rpt is not None and Path(h1_twas_rpt).exists() and os.path.getsize(h1_twas_rpt) > 0:
        twas_df = pd.read_table(h1_twas_rpt, usecols=[GENE_ID_COL_NAME, 'TWAS.P'])
        twas_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        twas_df.drop(labels=twas_df[twas_df['TWAS.P'] >= tool_thresholds['twas']].index, inplace=True)
        twas_df.drop(columns=[col for col in twas_df.columns if col != GENE_ID_COL_NAME], inplace=True)
    else:
        twas_df = pd.DataFrame(columns=[GENE_ID_COL_NAME])

    result_df = __merge_tool_gene_id(generated_file_path, coloc_df, smr_df, fastenloc_df, predixcan_df, ecaviar_df,
                                     twas_df)
    pre_mean = os.path.join(os.path.dirname(output_figure_path), f'pre_mean.tsv')
    result_df.to_csv(pre_mean, sep='\t', header=True, index=False, na_rep='NA')
    result_df.set_index(GENE_ID_COL_NAME, inplace=True)
    tools = [tool for tool, _, _ in TOOL_SIG_COL_INFO]
    count = []
    means = []
    stds = []
    for n in range(1, len(tools) + 1):
        count.append(n)
        tp_count = []
        for com in combinations(tools, n):
            uniq_df = pd.concat([result_df[tool] for tool in com])
            uniq_df.drop_duplicates(inplace=True)
            tp_count.append(uniq_df.notna().sum())
        tp_count_df = pd.DataFrame(tp_count, columns=['cnt'])
        means.append(tp_count_df['cnt'].mean())
        stds.append(0 if tp_count_df.shape[0] == 1 else tp_count_df['cnt'].std(ddof=0))
    print(f'count:{count}\nmeans:{means}\nstds:{stds}')
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    rects = ax.bar(count, means, 0.5, label='Mean')
    ax.bar_label(rects, fmt='{:0.1f}')
    ax.errorbar(count, means, yerr=stds, fmt=',', ecolor='black', capsize=5)
    ax.set_ylabel('Mean')
    ax.set_xlabel('Num of tools')
    ax.set_title(f'Mean of different combinations of tools(Total:{result_df.shape[0]})')
    # ax.set_xticks(x, count)
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    plt.close()


def __read_tool_info_align_to_prob(rpt, tool_name, sig_col_name, sig_type):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None or (not os.path.exists(rpt)) or os.path.getsize(rpt) <= 0:
        return None
    rpt_df = pd.read_table(rpt, usecols=[sig_col_name, GENE_ID_COL_NAME])
    rpt_df.sort_values(sig_col_name, ascending=sig_type == RESULT_TYPE_PVAL, inplace=True)
    rpt_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    if sig_type == RESULT_TYPE_PVAL:
        rpt_df[sig_col_name] = rpt_df[sig_col_name] * -1
    rpt_df.rename(columns={sig_col_name: tool_name}, inplace=True)
    return rpt_df


def plot_spearman_heatmap(rpts, output_figure_path=None, genetic_model='H1', tested_gene_df=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    df_list = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        tool_rpt = rpts.get(tool)
        if tool_rpt is None or len(tool_rpt) == 0:
            continue
        df_list.append(
            __read_tool_info_align_to_prob(tool_rpt, tool_name=tool, sig_col_name=sig_column, sig_type=sig_type))
    if tested_gene_df is None:
        # merge all genes
        gene_list = []
        for rpt_df in [df for df in df_list if df is not None and df.shape[0] > 0]:
            gene_list.append(rpt_df[[GENE_ID_COL_NAME]])
        gene_list_df = pd.concat(gene_list)
        del gene_list
        gene_list_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    else:
        gene_list_df = tested_gene_df[[GENE_ID_COL_NAME]]
    result_df = gene_list_df
    for rpt_df in [df for df in df_list if df is not None and df.shape[0] > 0]:
        result_df = pd.merge(left=result_df, right=rpt_df,
                             left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                             how='outer')
    # result_df[tool] is aligned to prob (for pval, result_df[tool] = pval * -1)
    # so fill NA with 0 for prob and -1 for pval
    for tool, _, sig_type in TOOL_SIG_COL_INFO:
        if tool in result_df.columns:
            result_df.loc[result_df[tool].isna(), tool] = 0 if sig_type == RESULT_TYPE_PROB else -1
    if result_df is None:
        return
    pre_heatmap = os.path.join(os.path.dirname(output_figure_path), f'{genetic_model}_pre_heatmap.tsv')
    result_df.to_csv(pre_heatmap, sep='\t', header=True, index=False, na_rep='NA')
    result_df.drop(columns=GENE_ID_COL_NAME, inplace=True)
    spearman_corr = result_df.corr(method='spearman')
    sns.set(font_scale=1.2)
    cluster_grid = sns.clustermap(spearman_corr, vmin=0, vmax=1, annot=True, cmap="vlag", fmt='.3f',
                                  cbar_kws=dict(ticks=[0, .2, .4, .6, .8, 1.0]))
    # heatmap.set_title('Spearman Correlation Heatmap')
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    plt.close()


def plot_upset(generated_file_path, h1_rpt_obj=None, sec_rpt_obj=None,
               sec_causal_type=1, output_figure_path=None, min_subset_size=3, sort_by='cardinality'):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    df_dict = {}
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        sec_rpt = sec_rpt_obj.get(tool)
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            tool_df = prepare_plot_data(generated_file_path, h1_rpt, sec_rpt, sec_causal_type,
                                        rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                        rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                        tool=tool)
            if tool_df.empty:
                continue
            df_dict[tool] = tool_df
    if len(df_dict) == 0:
        print('all report is empty, nothing to do')
        return
    threshold_param = {}
    ranking_df = None
    for tool, rpt_df in df_dict.items():
        rpt_df.drop(columns=[is_positive_col_name, prob_col_name], inplace=True)
        tool_pre_threshold_file = os.path.join(os.path.dirname(output_figure_path), f'{tool}_pre_threshold.tsv')
        rpt_df.to_csv(tool_pre_threshold_file, sep='\t', index=False, na_rep='NA')
        threshold_param[tool] = tool_pre_threshold_file
        if ranking_df is None:
            ranking_df = rpt_df
        else:
            ranking_df = pd.merge(left=ranking_df, right=rpt_df,
                                  left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                                  how='outer')
    tool_thresholds = calc.calc_threshold(rpt_obj=threshold_param, work_dir=os.path.dirname(output_figure_path))
    for _, tool_pre_threshold_file in threshold_param.items():
        os.remove(tool_pre_threshold_file)
    print(f'Thresholds: {tool_thresholds}')
    tool_cols = []
    tool_positive_cols = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        if sig_column not in ranking_df.columns:
            continue
        ranking_df.rename(columns={sig_column: tool}, inplace=True)
        # thresholds are covert to tool specific!
        if tool == 'smr':
            ranking_df[f'{tool}_positive'] = (ranking_df[tool] < tool_thresholds[tool]) & (
                    ranking_df['p_HEIDI'] > 0.05)
        elif sig_type == RESULT_TYPE_PROB:
            ranking_df[f'{tool}_positive'] = ranking_df[tool] > tool_thresholds[tool]
        elif sig_type == RESULT_TYPE_PVAL:
            ranking_df[f'{tool}_positive'] = ranking_df[tool] < tool_thresholds[tool]
        else:
            raise ValueError(f'Tool {tool} is not recognized, what is its threshold for it?')
        tool_cols.append(tool)
        tool_positive_cols.append(f'{tool}_positive')

    tool_val_cols = [f'{tool}_val' for tool in tool_cols]
    val_mapper = {k: v for k, v in zip(tool_cols, tool_val_cols)}
    ranking_df.rename(columns=val_mapper, inplace=True)
    positive_mapper = {k: v for k, v in zip(tool_positive_cols, tool_cols)}
    ranking_df.rename(columns=positive_mapper, inplace=True)
    pre_upset = os.path.join(os.path.dirname(output_figure_path), f'pre_upset.tsv')
    ranking_df.to_csv(pre_upset, sep='\t', index=False, na_rep='NA')
    plt.figure().clear()
    plot(from_indicators(tool_cols, data=ranking_df[tool_cols]),
         show_counts=True,
         min_subset_size=min_subset_size,
         sort_by=sort_by)
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--generated_file', dest='generated_file',
                        required=True,
                        help='Simulation generated file list file')
    parser.add_argument('--h1_coloc_rpt', dest='h1_coloc_rpt',
                        help='H1 coloc result file')
    parser.add_argument('--h1_smr_rpt', dest='h1_smr_rpt',
                        help='H1 smr result file')
    # parser.add_argument('--h1_jlim_rpt', dest='h1_jlim_rpt',
    #                     help='H1 jlim result file')
    parser.add_argument('--h1_fastenloc_rpt', dest='h1_fastenloc_rpt',
                        help='H1 fastenloc result file')
    parser.add_argument('--h1_predixcan_rpt', dest='h1_predixcan_rpt',
                        help='H1 predixcan result file')
    parser.add_argument('--h1_ecaviar_rpt', dest='h1_ecaviar_rpt',
                        help='H1 ecaviar result file')
    parser.add_argument('--h1_twas_rpt', dest='h1_twas_rpt',
                        help='H1 twas result file')

    parser.add_argument('--sec_coloc_rpt', dest='sec_coloc_rpt',
                        help='H0/H2 coloc result file')
    parser.add_argument('--sec_smr_rpt', dest='sec_smr_rpt',
                        help='H0/H2 smr result file')
    # parser.add_argument('--sec_jlim_rpt', dest='sec_jlim_rpt',
    #                     help='H0/H2 jlim result file')
    parser.add_argument('--sec_fastenloc_rpt', dest='sec_fastenloc_rpt',
                        help='H0/H2 fastenloc result file')
    parser.add_argument('--sec_predixcan_rpt', dest='sec_predixcan_rpt',
                        help='H0/H2 predixcan result file')
    parser.add_argument('--sec_ecaviar_rpt', dest='sec_ecaviar_rpt',
                        help='H0/H2 ecaviar result file')
    parser.add_argument('--sec_twas_rpt', dest='sec_twas_rpt',
                        help='H0/H2 twas result file')

    parser.add_argument('--sec_causal_type', dest='sec_causal_type', type=int, required=True, choices=range(0, 8),
                        help='H0/H2 causal type: 0:H0; 1:H1; 2:H2 r2<=0.4; 3:H2 0.4<r2<=0.7; 4:H2 0.7<r2<=0.9')
    parser.add_argument('--roc_output_figure_path', dest='roc_output_figure_path',
                        help='ROC figure output file path')
    parser.add_argument('--prc_output_figure_path', dest='prc_output_figure_path',
                        help='PRC figure output file path')
    args = parser.parse_args()
    print(f'Accepted args:\n {args}')
    plot_all_roc(args.generated_file,
                 {'coloc': args.h1_coloc_rpt,
                  'smr': args.h1_smr_rpt,
                  'fastenloc': args.h1_fastenloc_rpt,
                  'predixcan': args.h1_predixcan_rpt,
                  'ecaviar': args.h1_ecaviar_rpt,
                  'twas': args.h1_twas_rpt},
                 {'coloc': args.sec_coloc_rpt,
                  'smr': args.sec_smr_rpt,
                  'fastenloc': args.sec_fastenloc_rpt,
                  'predixcan': args.sec_predixcan_rpt,
                  'ecaviar': args.sec_ecaviar_rpt,
                  'twas': args.sec_twas_rpt},
                 args.sec_causal_type, args.roc_output_figure_path)
    plot_all_prc(args.generated_file,
                 {'coloc': args.h1_coloc_rpt,
                  'smr': args.h1_smr_rpt,
                  'fastenloc': args.h1_fastenloc_rpt,
                  'predixcan': args.h1_predixcan_rpt,
                  'ecaviar': args.h1_ecaviar_rpt,
                  'twas': args.h1_twas_rpt},
                 {'coloc': args.sec_coloc_rpt,
                  'smr': args.sec_smr_rpt,
                  'fastenloc': args.sec_fastenloc_rpt,
                  'predixcan': args.sec_predixcan_rpt,
                  'ecaviar': args.sec_ecaviar_rpt,
                  'twas': args.sec_twas_rpt},
                 args.sec_causal_type, args.prc_output_figure_path)
