from pathlib import Path
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from sklearn.metrics import roc_auc_score, roc_curve, average_precision_score, precision_recall_curve
import ranking.birra as birra
import ranking.rra as rra
import ranking.intact as intact

prob_col_name = 'PROB'
is_positive_col_name = 'IS_POSITIVE'


def retrieve_std_df(generated_list, sec_causal_type):
    generated_list_df = pd.read_table(generated_list,
                                      usecols=['gene', 'gwas_causal_snp', 'gwas_max_assoc_p', 'gwas_r2',
                                               'eqtl_max_assoc_p_1', f'eqtl_max_assoc_p_{sec_causal_type}'])
    positive_df = generated_list_df[(generated_list_df['gwas_max_assoc_p'] <= 1.0E-5) & (
            generated_list_df['gwas_r2'] >= 0.8) & (generated_list_df['eqtl_max_assoc_p_1'] <= 0.01)].copy()
    positive_df[is_positive_col_name] = 1
    positive_df['gene_id'] = positive_df['gene'] + '_1'

    negative_df = generated_list_df[(generated_list_df['gwas_max_assoc_p'] <= 1.0E-5) & (
            generated_list_df['gwas_r2'] >= 0.8) & (generated_list_df[
                                                        f'eqtl_max_assoc_p_{sec_causal_type}'] <= 0.01)].copy()
    negative_df[is_positive_col_name] = 0
    negative_df['gene_id'] = negative_df['gene'] + f'_{sec_causal_type}'

    std_df = pd.concat([positive_df[['gene_id', is_positive_col_name]], negative_df[['gene_id', is_positive_col_name]]])
    return std_df


def prepare_plot_data(generated_list, h1_report, sec_report, sec_causal_type,
                      rpt_prob_col_name=None, rpt_pval_col_name=None, tool=None):
    if rpt_prob_col_name is None and rpt_pval_col_name is None:
        raise ValueError('p-value column name and probability column name can not be both null')
    std_df = retrieve_std_df(generated_list, sec_causal_type)

    if rpt_prob_col_name is None:
        reading_cols = [rpt_pval_col_name, 'gene_id']
        if tool == 'predixcan':
            reading_cols.append('zscore')
        elif tool == 'twas':
            reading_cols.append('TWAS.Z')
        elif tool == 'smr':
            reading_cols.append('b_SMR')
            reading_cols.append('se_SMR')
        h1_report_df = pd.read_table(h1_report, usecols=reading_cols)
        h1_report_df.drop_duplicates(subset='gene_id', inplace=True)
        h1_report_df[prob_col_name] = 1 - h1_report_df[rpt_pval_col_name]
        if sec_report is not None and Path(sec_report).exists() and os.path.getsize(sec_report) > 0:
            sec_report_df = pd.read_table(sec_report, usecols=reading_cols)
            sec_report_df.drop_duplicates(subset='gene_id', inplace=True)
            # Convert p-value to probability. TODO this method is not good
            sec_report_df[prob_col_name] = 1 - sec_report_df[rpt_pval_col_name]
        else:
            sec_report_df = pd.DataFrame(columns=[rpt_pval_col_name, prob_col_name, 'gene_id'])
    else:
        h1_report_df = pd.read_table(h1_report, usecols=[rpt_prob_col_name, 'gene_id'])
        h1_report_df.drop_duplicates(subset='gene_id', inplace=True)
        h1_report_df[prob_col_name] = h1_report_df[rpt_prob_col_name]
        if sec_report is not None and Path(sec_report).exists() and os.path.getsize(sec_report) > 0:
            sec_report_df = pd.read_table(sec_report, usecols=[rpt_prob_col_name, 'gene_id'])
            sec_report_df.drop_duplicates(subset='gene_id', inplace=True)
            sec_report_df[prob_col_name] = sec_report_df[rpt_prob_col_name]
        else:
            sec_report_df = pd.DataFrame(columns=[rpt_prob_col_name, prob_col_name, 'gene_id'])
    h1_report_df['gene_id'] = h1_report_df['gene_id'] + '_1'
    sec_report_df['gene_id'] = sec_report_df['gene_id'] + f'_{sec_causal_type}'
    report_df = pd.concat([h1_report_df, sec_report_df])
    result_df = pd.merge(left=std_df, right=report_df,
                         left_on='gene_id', right_on='gene_id',
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
    result_df = prepare_plot_data(generated_list, h1_report, sec_report, sec_causal_type, rpt_prob_col_name,
                                  rpt_pval_col_name, tool=tool)
    plot_roc_curve(result_df[is_positive_col_name], result_df[prob_col_name], tool=tool)


def plot_roc_curve(true_y, y_prob, tool):
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    auc = roc_auc_score(true_y, y_prob)
    if tool is None:
        lbl = 'AUC = {:0.4f}'.format(auc)
    else:
        lbl = tool + ': AUC = {:0.4f}'
        lbl = lbl.format(auc)
    plt.plot(fpr, tpr, label=lbl)
    plt.legend(loc='lower right')


def plot_all_roc(generated_file_path,
                 h1_coloc_rpt, sec_coloc_rpt,
                 h1_smr_rpt, sec_smr_rpt,
                 h1_jlim_rpt, sec_jlim_rpt,
                 h1_fastenloc_rpt=None, sec_fastenloc_rpt=None,
                 h1_predixcan_rpt=None, sec_predixcan_rpt=None,
                 h1_ecaviar_rpt=None, sec_ecaviar_rpt=None,
                 h1_twas_rpt=None, sec_twas_rpt=None,
                 sec_causal_type=1,
                 output_figure_path=None):
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
    if h1_coloc_rpt is not None and Path(h1_coloc_rpt).exists() and os.path.getsize(h1_coloc_rpt) > 0:
        plot_single_roc(generated_file_path, h1_coloc_rpt, sec_coloc_rpt, sec_causal_type,
                        rpt_prob_col_name='overall_H4', rpt_pval_col_name=None,
                        tool='coloc')
    if h1_smr_rpt is not None and Path(h1_smr_rpt).exists() and os.path.getsize(h1_smr_rpt) > 0:
        plot_single_roc(generated_file_path, h1_smr_rpt, sec_smr_rpt, sec_causal_type,
                        rpt_prob_col_name=None, rpt_pval_col_name='p_SMR', tool='smr')
    if h1_jlim_rpt is not None and Path(h1_jlim_rpt).exists() and os.path.getsize(h1_jlim_rpt) > 0:
        plot_single_roc(generated_file_path, h1_jlim_rpt, sec_jlim_rpt, sec_causal_type,
                        rpt_prob_col_name=None, rpt_pval_col_name='pvalue', tool='jlim')
    if h1_fastenloc_rpt is not None and Path(h1_fastenloc_rpt).exists() and os.path.getsize(h1_fastenloc_rpt) > 0:
        plot_single_roc(generated_file_path, h1_fastenloc_rpt, sec_fastenloc_rpt, sec_causal_type,
                        rpt_prob_col_name='RCP', rpt_pval_col_name=None,
                        tool='fastenloc')
    if h1_predixcan_rpt is not None and Path(h1_predixcan_rpt).exists() and os.path.getsize(h1_predixcan_rpt) > 0:
        plot_single_roc(generated_file_path, h1_predixcan_rpt, sec_predixcan_rpt, sec_causal_type,
                        rpt_prob_col_name=None, rpt_pval_col_name='pvalue', tool='predixcan')
    if h1_ecaviar_rpt is not None and Path(h1_ecaviar_rpt).exists() and os.path.getsize(h1_ecaviar_rpt) > 0:
        plot_single_roc(generated_file_path, h1_ecaviar_rpt, sec_ecaviar_rpt, sec_causal_type,
                        rpt_prob_col_name='clpp', rpt_pval_col_name=None,
                        tool='ecaviar')
    if h1_twas_rpt is not None and Path(h1_twas_rpt).exists() and os.path.getsize(h1_twas_rpt) > 0:
        plot_single_roc(generated_file_path, h1_twas_rpt, sec_twas_rpt, sec_causal_type,
                        rpt_prob_col_name=None, rpt_pval_col_name='TWAS.P',
                        tool='twas')
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()


def plot_single_prc(generated_list, h1_report, sec_report, sec_causal_type,
                    rpt_prob_col_name=None, rpt_pval_col_name=None, tool=None):
    result_df = prepare_plot_data(generated_list, h1_report, sec_report, sec_causal_type, rpt_prob_col_name,
                                  rpt_pval_col_name)
    plot_prc_curve(result_df[is_positive_col_name], result_df[prob_col_name], tool=tool)


def plot_all_prc(generated_file_path,
                 h1_coloc_rpt, sec_coloc_rpt,
                 h1_smr_rpt, sec_smr_rpt,
                 h1_jlim_rpt, sec_jlim_rpt,
                 h1_fastenloc_rpt=None, sec_fastenloc_rpt=None,
                 h1_predixcan_rpt=None, sec_predixcan_rpt=None,
                 h1_ecaviar_rpt=None, sec_ecaviar_rpt=None,
                 sec_causal_type=1,
                 output_figure_path=None):
    plt.figure().clear()
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    if h1_coloc_rpt is not None and Path(h1_coloc_rpt).exists() and os.path.getsize(h1_coloc_rpt) > 0:
        plot_single_prc(generated_file_path, h1_coloc_rpt, sec_coloc_rpt, sec_causal_type,
                        rpt_prob_col_name='overall_H4', rpt_pval_col_name=None, tool='coloc')
    if h1_smr_rpt is not None and Path(h1_smr_rpt).exists() and os.path.getsize(h1_smr_rpt) > 0:
        plot_single_prc(generated_file_path, h1_smr_rpt, sec_smr_rpt, sec_causal_type,
                        rpt_prob_col_name=None, rpt_pval_col_name='p_SMR', tool='smr')
    if h1_jlim_rpt is not None and Path(h1_jlim_rpt).exists() and os.path.getsize(h1_jlim_rpt) > 0:
        plot_single_prc(generated_file_path, h1_jlim_rpt, sec_jlim_rpt, sec_causal_type,
                        rpt_prob_col_name=None, rpt_pval_col_name='pvalue', tool='jlim')
    if h1_fastenloc_rpt is not None and Path(h1_fastenloc_rpt).exists() and os.path.getsize(h1_fastenloc_rpt) > 0:
        plot_single_prc(generated_file_path, h1_fastenloc_rpt, sec_fastenloc_rpt, sec_causal_type,
                        rpt_prob_col_name='RCP', rpt_pval_col_name=None, tool='fastenloc')
    if h1_predixcan_rpt is not None and Path(h1_predixcan_rpt).exists() and os.path.getsize(h1_predixcan_rpt) > 0:
        plot_single_prc(generated_file_path, h1_predixcan_rpt, sec_predixcan_rpt, sec_causal_type,
                        rpt_prob_col_name=None, rpt_pval_col_name='pvalue', tool='predixcan')
    if h1_ecaviar_rpt is not None and Path(h1_ecaviar_rpt).exists() and os.path.getsize(h1_ecaviar_rpt) > 0:
        plot_single_prc(generated_file_path, h1_ecaviar_rpt, sec_ecaviar_rpt, sec_causal_type,
                        rpt_prob_col_name='clpp', rpt_pval_col_name=None, tool='eCaviar')
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()


def plot_prc_curve(true_y, y_prob, tool):
    precision, recall, thresholds = precision_recall_curve(true_y, y_prob)
    ap = average_precision_score(true_y, y_prob)
    if tool is None:
        lbl = 'AP = {:0.4f}'.format(ap)
    else:
        lbl = tool + ': AP = {:0.4f}'
        lbl = lbl.format(ap)
    plt.plot(recall, precision, label=lbl)
    plt.legend(loc='upper right')


def plot_all_against_ensemble_roc(generated_file_path,
                                  h1_coloc_rpt, sec_coloc_rpt,
                                  h1_smr_rpt, sec_smr_rpt,
                                  h1_jlim_rpt, sec_jlim_rpt,
                                  h1_fastenloc_rpt=None, sec_fastenloc_rpt=None,
                                  h1_predixcan_rpt=None, sec_predixcan_rpt=None,
                                  h1_ecaviar_rpt=None, sec_ecaviar_rpt=None,
                                  h1_twas_rpt=None, sec_twas_rpt=None,
                                  sec_causal_type=1,
                                  output_figure_path=None, output_dir=''):
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
    if h1_coloc_rpt is not None and Path(h1_coloc_rpt).exists() and os.path.getsize(h1_coloc_rpt) > 0:
        coloc_df = prepare_plot_data(generated_file_path, h1_coloc_rpt, sec_coloc_rpt, sec_causal_type,
                                     rpt_prob_col_name='overall_H4', rpt_pval_col_name=None, tool='coloc')
        plot_roc_curve(coloc_df[is_positive_col_name], coloc_df[prob_col_name], tool='coloc')
    else:
        coloc_df = None
    if h1_smr_rpt is not None and Path(h1_smr_rpt).exists() and os.path.getsize(h1_smr_rpt) > 0:
        smr_df = prepare_plot_data(generated_file_path, h1_smr_rpt, sec_smr_rpt, sec_causal_type,
                                   rpt_prob_col_name=None, rpt_pval_col_name='p_SMR', tool='smr')
        plot_roc_curve(smr_df[is_positive_col_name], smr_df[prob_col_name], tool='SMR')
    else:
        smr_df = None
    if h1_jlim_rpt is not None and Path(h1_jlim_rpt).exists() and os.path.getsize(h1_jlim_rpt) > 0:
        jlim_df = prepare_plot_data(generated_file_path, h1_jlim_rpt, sec_jlim_rpt, sec_causal_type,
                                    rpt_prob_col_name=None, rpt_pval_col_name='pvalue', tool='jlim')
        plot_roc_curve(jlim_df[is_positive_col_name], jlim_df[prob_col_name], tool='JLIM')
    else:
        jlim_df = None
    if h1_fastenloc_rpt is not None and Path(h1_fastenloc_rpt).exists() and os.path.getsize(h1_fastenloc_rpt) > 0:
        fastenloc_df = prepare_plot_data(generated_file_path, h1_fastenloc_rpt, sec_fastenloc_rpt, sec_causal_type,
                                         rpt_prob_col_name='RCP', rpt_pval_col_name=None, tool='fastenloc')
        plot_roc_curve(fastenloc_df[is_positive_col_name], fastenloc_df[prob_col_name], tool='fastenloc')
    else:
        fastenloc_df = None

    if h1_predixcan_rpt is not None and Path(h1_predixcan_rpt).exists() and os.path.getsize(h1_predixcan_rpt) > 0:
        predixcan_df = prepare_plot_data(generated_file_path, h1_predixcan_rpt, sec_predixcan_rpt, sec_causal_type,
                                         rpt_prob_col_name=None, rpt_pval_col_name='pvalue', tool='predixcan')
        plot_roc_curve(predixcan_df[is_positive_col_name], predixcan_df[prob_col_name], tool='predixcan')
    else:
        predixcan_df = None

    if h1_ecaviar_rpt is not None and Path(h1_ecaviar_rpt).exists() and os.path.getsize(h1_ecaviar_rpt) > 0:
        ecaviar_df = prepare_plot_data(generated_file_path, h1_ecaviar_rpt, sec_ecaviar_rpt, sec_causal_type,
                                       rpt_prob_col_name='clpp', rpt_pval_col_name=None, tool='ecaviar')
        plot_roc_curve(ecaviar_df[is_positive_col_name], ecaviar_df[prob_col_name], tool='eCaviar')
    else:
        ecaviar_df = None

    if h1_twas_rpt is not None and Path(h1_twas_rpt).exists() and os.path.getsize(h1_twas_rpt) > 0:
        twas_df = prepare_plot_data(generated_file_path, h1_twas_rpt, sec_twas_rpt, sec_causal_type,
                                    rpt_prob_col_name=None, rpt_pval_col_name='TWAS.P', tool='twas')
        plot_roc_curve(twas_df[is_positive_col_name], twas_df[prob_col_name], tool='TWAS')
    else:
        twas_df = None

    if coloc_df is None and smr_df is None and jlim_df is None and fastenloc_df is None and \
            predixcan_df is None and ecaviar_df is None and twas_df is None:
        print('all report is empty, nothing to do')
        return

    if coloc_df is not None:
        coloc_ranking_file = os.path.join(output_dir, 'coloc_result.tsv')
        coloc_df.to_csv(coloc_ranking_file, sep='\t', header=True, index=False)
    else:
        coloc_ranking_file = None
    if smr_df is not None:
        smr_ranking_file = os.path.join(output_dir, 'smr_result.tsv')
        smr_df.to_csv(smr_ranking_file, sep='\t', header=True, index=False)
    else:
        smr_ranking_file = None
    if jlim_df is not None:
        jlim_ranking_file = os.path.join(output_dir, 'jlim_result.tsv')
        jlim_df.to_csv(jlim_ranking_file, sep='\t', header=True, index=False)
    else:
        jlim_ranking_file = None
    if fastenloc_df is not None:
        fastenloc_ranking_file = os.path.join(output_dir, 'fastenloc_result.tsv')
        fastenloc_df.to_csv(fastenloc_ranking_file, sep='\t', header=True, index=False)
    else:
        fastenloc_ranking_file = None

    if predixcan_df is not None:
        predixcan_ranking_file = os.path.join(output_dir, 'predixcan_result.tsv')
        predixcan_df.to_csv(predixcan_ranking_file, sep='\t', header=True, index=False)
    else:
        predixcan_ranking_file = None

    if ecaviar_df is not None:
        ecaviar_ranking_file = os.path.join(output_dir, 'ecaviar_result.tsv')
        ecaviar_df.to_csv(ecaviar_ranking_file, sep='\t', header=True, index=False)
    else:
        ecaviar_ranking_file = None

    if twas_df is not None:
        twas_ranking_file = os.path.join(output_dir, 'twas_result.tsv')
        twas_df.to_csv(twas_ranking_file, sep='\t', header=True, index=False)
    else:
        twas_ranking_file = None

    ensemble_output_file = os.path.join(output_dir, 'birra.tsv')
    rpt_obj = {
        'coloc': coloc_ranking_file,
        'smr': smr_ranking_file,
        'jlim': jlim_ranking_file,
        'fastenloc': fastenloc_ranking_file,
        'predixcan': predixcan_ranking_file,
        'ecaviar': ecaviar_ranking_file,
        'twas': twas_ranking_file
    }
    birra.run_ranking(output_file_path=ensemble_output_file, rpt=rpt_obj)
    ranking_result_df = pd.read_table(ensemble_output_file, usecols=['gene_id', 'birra_ranking'])
    # Convert birra ranking to probability. TODO this method is not good
    # ranking_result_df[prob_col_name] = 1 - ranking_result_df['birra_ranking'] / ranking_result_df['birra_ranking'].max()
    std_df = retrieve_std_df(generated_file_path, sec_causal_type)
    # ranking_result_df = pd.merge(left=ranking_result_df, right=std_df,
    #                              left_on='gene_id', right_on='gene_id',
    #                              how='left')
    # tp_na_bool_series = ranking_result_df[is_positive_col_name].isna()
    # ranking_result_df.loc[ranking_result_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(ranking_result_df[is_positive_col_name], ranking_result_df[prob_col_name], tool='birra_ranking')
    # --------
    rra_geo_output_file = os.path.join(output_dir, 'rra_geo.tsv')
    rra.run_ranking(output_file_path=rra_geo_output_file, rpt=rpt_obj, sample_size=191647, method='GEO')
    rgeo_df = pd.read_table(rra_geo_output_file, usecols=['gene_id', 'geo_p_value'])
    # Convert rra p_value to probability. TODO this method is not good
    rgeo_df[prob_col_name] = 1 - rgeo_df['geo_p_value']
    rgeo_df = pd.merge(left=rgeo_df, right=std_df,
                       left_on='gene_id', right_on='gene_id',
                       how='left')
    tp_na_bool_series = rgeo_df[is_positive_col_name].isna()
    rgeo_df.loc[rgeo_df[tp_na_bool_series].index, is_positive_col_name] = 0
    plot_roc_curve(rgeo_df[is_positive_col_name], rgeo_df[prob_col_name], tool='rGEO')
    # --------
    # rra_stuart_output_file = os.path.join(output_dir, 'rra_stuart.tsv')
    # rra.run_ranking(output_file_path=rra_stuart_output_file, rpt=rpt_obj, sample_size=191647, method='stuart')
    # stuart_df = pd.read_table(rra_stuart_output_file, usecols=['gene_id', 'stuart_p_value'])
    # # Convert rra p_value to probability. TODO this method is not good
    # stuart_df[prob_col_name] = 1 - stuart_df['stuart_p_value']
    # stuart_df = pd.merge(left=stuart_df, right=std_df,
    #                      left_on='gene_id', right_on='gene_id',
    #                      how='left')
    # tp_na_bool_series = stuart_df[is_positive_col_name].isna()
    # stuart_df.loc[stuart_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(stuart_df[is_positive_col_name], stuart_df[prob_col_name], tool='stuart')
    # --------xgboost lr
    # train_output_file = os.path.join(output_dir, 'train_xgb_lr.tsv')
    # machine_learning.run_ranking(rpt=rpt_obj, output_file_path=train_output_file)
    # train_df = pd.read_table(train_output_file, usecols=['gene_id', 'ml_probability'])
    # train_df = pd.merge(left=train_df, right=std_df,
    #                     left_on='gene_id', right_on='gene_id',
    #                     how='left')
    # tp_na_bool_series = train_df[is_positive_col_name].isna()
    # train_df.loc[train_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_df[is_positive_col_name], train_df['ml_probability'], tool='ensemble')
    # --------xgboost rf
    # train_xgb_rf_output_file = os.path.join(output_dir, 'train_xgb_rf.tsv')
    # machine_learning.run_ranking_xgb_rf(rpt=rpt_obj, output_file_path=train_xgb_rf_output_file)
    # train_xgb_rf_df = pd.read_table(train_xgb_rf_output_file, usecols=['gene_id', 'ml_probability'])
    # train_xgb_rf_df = pd.merge(left=train_xgb_rf_df, right=std_df,
    #                            left_on='gene_id', right_on='gene_id',
    #                            how='left')
    # tp_na_bool_series = train_xgb_rf_df[is_positive_col_name].isna()
    # train_xgb_rf_df.loc[train_xgb_rf_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_xgb_rf_df[is_positive_col_name], train_xgb_rf_df['ml_probability'], tool='ensemble xgb rf')
    # --------scikit-learn svm
    # train_svm_output_file = os.path.join(output_dir, 'train_svm.tsv')
    # machine_learning.run_ranking_svm(rpt=rpt_obj, output_file_path=train_svm_output_file)
    # train_svm_df = pd.read_table(train_svm_output_file, usecols=['gene_id', 'ml_probability'])
    # train_svm_df = pd.merge(left=train_svm_df, right=std_df,
    #                         left_on='gene_id', right_on='gene_id',
    #                         how='left')
    # tp_na_bool_series = train_svm_df[is_positive_col_name].isna()
    # train_svm_df.loc[train_svm_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_svm_df[is_positive_col_name], train_svm_df['ml_probability'], tool='ensemble svm')
    # --------scikit-learn random forest
    # train_rf_output_file = os.path.join(output_dir, 'train_rf.tsv')
    # machine_learning.run_ranking_rf(rpt=rpt_obj, output_file_path=train_rf_output_file)
    # train_rf_df = pd.read_table(train_rf_output_file, usecols=['gene_id', 'ml_probability'])
    # train_rf_df = pd.merge(left=train_rf_df, right=std_df,
    #                        left_on='gene_id', right_on='gene_id',
    #                        how='left')
    # tp_na_bool_series = train_rf_df[is_positive_col_name].isna()
    # train_rf_df.loc[train_rf_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_rf_df[is_positive_col_name], train_rf_df['ml_probability'], tool='ensemble rf')
    # --------scikit-learn logistic regression
    # train_lr_output_file = os.path.join(output_dir, 'train_lr.tsv')
    # machine_learning.run_ranking_lr(rpt=rpt_obj, output_file_path=train_lr_output_file)
    # train_lr_df = pd.read_table(train_lr_output_file, usecols=['gene_id', 'ml_probability'])
    # train_lr_df = pd.merge(left=train_lr_df, right=std_df,
    #                        left_on='gene_id', right_on='gene_id',
    #                        how='left')
    # tp_na_bool_series = train_lr_df[is_positive_col_name].isna()
    # train_lr_df.loc[train_lr_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(train_lr_df[is_positive_col_name], train_lr_df['ml_probability'], tool='ensemble lr')
    # --------intact
    # intact_output_file = os.path.join(output_dir, 'intact_linear.tsv')
    # intact.run_ranking(rpt=rpt_obj, output_file_path=intact_output_file, prior_fun='linear')
    # intact_df = pd.read_table(intact_output_file, usecols=['gene_id', 'intact_probability'])
    # intact_df = pd.merge(left=intact_df, right=std_df,
    #                      left_on='gene_id', right_on='gene_id',
    #                      how='left')
    # tp_na_bool_series = intact_df[is_positive_col_name].isna()
    # intact_df.loc[intact_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # intact_df.loc[intact_df['intact_probability'].isna(), 'intact_probability'] = 0
    # plot_roc_curve(intact_df[is_positive_col_name], intact_df['intact_probability'], tool='intact linear')
    # --------intact_hybrid
    intact_hybrid_output_file = os.path.join(output_dir, 'intact_hybrid.tsv')
    intact.run_ranking(rpt=rpt_obj, output_file_path=intact_hybrid_output_file, prior_fun='hybrid')
    intact_hybrid_df = pd.read_table(intact_hybrid_output_file, usecols=['gene_id', 'intact_probability'])
    intact_hybrid_df = pd.merge(left=intact_hybrid_df, right=std_df,
                                left_on='gene_id', right_on='gene_id',
                                how='left')
    tp_na_bool_series = intact_hybrid_df[is_positive_col_name].isna()
    intact_hybrid_df.loc[intact_hybrid_df[tp_na_bool_series].index, is_positive_col_name] = 0
    intact_hybrid_df.loc[intact_hybrid_df['intact_probability'].isna(), 'intact_probability'] = 0
    plot_roc_curve(intact_hybrid_df[is_positive_col_name], intact_hybrid_df['intact_probability'], tool='intact')
    # --------intact_step
    # intact_step_output_file = os.path.join(output_dir, 'intact_step.tsv')
    # intact.run_ranking(rpt=rpt_obj, output_file_path=intact_step_output_file, prior_fun='step')
    # intact_step_df = pd.read_table(intact_step_output_file, usecols=['gene_id', 'intact_probability'])
    # intact_step_df = pd.merge(left=intact_step_df, right=std_df,
    #                           left_on='gene_id', right_on='gene_id',
    #                           how='left')
    # tp_na_bool_series = intact_step_df[is_positive_col_name].isna()
    # intact_step_df.loc[intact_step_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # intact_step_df.loc[intact_step_df['intact_probability'].isna(), 'intact_probability'] = 0
    # plot_roc_curve(intact_step_df[is_positive_col_name], intact_step_df['intact_probability'], tool='intact step')
    # --------intact_expit
    # intact_expit_output_file = os.path.join(output_dir, 'intact_expit.tsv')
    # intact.run_ranking(rpt=rpt_obj, output_file_path=intact_expit_output_file, prior_fun='expit')
    # intact_expit_df = pd.read_table(intact_expit_output_file, usecols=['gene_id', 'intact_probability'])
    # intact_expit_df = pd.merge(left=intact_expit_df, right=std_df,
    #                            left_on='gene_id', right_on='gene_id',
    #                            how='left')
    # tp_na_bool_series = intact_expit_df[is_positive_col_name].isna()
    # intact_expit_df.loc[intact_expit_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # intact_expit_df.loc[intact_expit_df['intact_probability'].isna(), 'intact_probability'] = 0
    # plot_roc_curve(intact_expit_df[is_positive_col_name], intact_expit_df['intact_probability'], tool='intact expit')
    # ---------simple ranking
    # simple_ranking_result_df = pd.read_table(intact_expit_output_file, usecols=['gene_id', 'avg_ranking'])
    # # Convert simple mean ranking to probability. TODO this method is not good
    # simple_ranking_result_df[prob_col_name] = 1 - simple_ranking_result_df['avg_ranking'] / simple_ranking_result_df[
    #     'avg_ranking'].max()
    # std_df = retrieve_std_df(generated_file_path, sec_causal_type)
    # simple_ranking_result_df = pd.merge(left=simple_ranking_result_df, right=std_df,
    #                                     left_on='gene_id', right_on='gene_id',
    #                                     how='left')
    # tp_na_bool_series = simple_ranking_result_df[is_positive_col_name].isna()
    # simple_ranking_result_df.loc[simple_ranking_result_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(simple_ranking_result_df[is_positive_col_name], simple_ranking_result_df[prob_col_name],
    #                tool='avg ranking')
    # --------
    # --------pu learning
    # pu_learning_output_file = os.path.join(output_dir, 'train_pu.tsv')
    # machine_learning.run_ranking_pu(rpt=rpt_obj, output_file_path=pu_learning_output_file)
    # pu_learning_df = pd.read_table(pu_learning_output_file, usecols=['gene_id', 'ml_probability'])
    # pu_learning_df = pd.merge(left=pu_learning_df, right=std_df,
    #                           left_on='gene_id', right_on='gene_id',
    #                           how='left')
    # tp_na_bool_series = pu_learning_df[is_positive_col_name].isna()
    # pu_learning_df.loc[pu_learning_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_roc_curve(pu_learning_df[is_positive_col_name], pu_learning_df['ml_probability'], tool='pu learning')
    # --------
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--generated_file', dest='generated_file',
                        required=True,
                        help='Simulation generated file list file')
    parser.add_argument('--h1_coloc_rpt', dest='h1_coloc_rpt',
                        help='H1 coloc result file')
    parser.add_argument('--h1_smr_rpt', dest='h1_smr_rpt',
                        help='H1 smr result file')
    parser.add_argument('--h1_jlim_rpt', dest='h1_jlim_rpt',
                        help='H1 jlim result file')
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
    parser.add_argument('--sec_jlim_rpt', dest='sec_jlim_rpt',
                        help='H0/H2 jlim result file')
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
    plot_all_roc(args.generated_file, args.h1_coloc_rpt, args.sec_coloc_rpt, args.h1_smr_rpt, args.sec_smr_rpt,
                 args.h1_jlim_rpt, args.sec_jlim_rpt, args.h1_fastenloc_rpt, args.sec_fastenloc_rpt,
                 args.h1_predixcan_rpt, args.sec_predixcan_rpt, args.h1_ecaviar_rpt, args.sec_ecaviar_rpt,
                 args.sec_causal_type, args.roc_output_figure_path)
    plot_all_prc(args.generated_file, args.h1_coloc_rpt, args.sec_coloc_rpt, args.h1_smr_rpt, args.sec_smr_rpt,
                 args.h1_jlim_rpt, args.sec_jlim_rpt, args.h1_fastenloc_rpt, args.sec_fastenloc_rpt,
                 args.h1_predixcan_rpt, args.sec_predixcan_rpt, args.h1_ecaviar_rpt, args.sec_ecaviar_rpt,
                 args.sec_causal_type, args.prc_output_figure_path)

    # ============testing============
    # _output_dir = '/Users/admin/Downloads/coloc_results'
    #
    # _generated_file = '/Users/admin/Downloads/coloc_results/generated_20221102031352.tsv'
    # _h1_coloc_rpt = '/Users/admin/Downloads/coloc_results/h1/coloc_output_20221212140605.tsv.gz'
    # if _h1_coloc_rpt is not None and (_h1_coloc_rpt == '' or _h1_coloc_rpt.lower() == "none"):
    #     _h1_coloc_rpt = None
    #
    # _h1_smr_rpt = '/Users/admin/Downloads/coloc_results/h1/smr_output_20221212141459.tsv.gz'
    # if _h1_smr_rpt is not None and (_h1_smr_rpt == '' or _h1_smr_rpt.lower() == "none"):
    #     _h1_smr_rpt = None
    #
    # _h1_jlim_rpt = None
    # if _h1_jlim_rpt is not None and (_h1_jlim_rpt == '' or _h1_jlim_rpt.lower() == "none"):
    #     _h1_jlim_rpt = None
    #
    # _h1_fastenloc_rpt = '/Users/admin/Downloads/coloc_results/h1/fastenloc_output_20221212153042.tsv.gz'
    # if _h1_fastenloc_rpt is not None and (_h1_fastenloc_rpt == '' or _h1_fastenloc_rpt.lower() == "none"):
    #     _h1_fastenloc_rpt = None
    #
    # _h1_predixcan_rpt = '/Users/admin/Downloads/coloc_results/h1/predixcan_output_20221212140610.tsv.gz'
    # if _h1_predixcan_rpt is not None and (_h1_predixcan_rpt == '' or _h1_predixcan_rpt.lower() == "none"):
    #     _h1_predixcan_rpt = None
    #
    # _h1_ecaviar_rpt = '/Users/admin/Downloads/coloc_results/h1/ecaviar_output_20221212141418.tsv.gz'
    # if _h1_ecaviar_rpt is not None and (_h1_ecaviar_rpt == '' or _h1_ecaviar_rpt.lower() == "none"):
    #     _h1_ecaviar_rpt = None
    #
    # _h1_twas_rpt = '/Users/admin/Downloads/coloc_results/h1/twas_output_20221223222006.tsv.gz'
    # if _h1_twas_rpt is not None and (_h1_twas_rpt == '' or _h1_twas_rpt.lower() == "none"):
    #     _h1_twas_rpt = None
    # ----------------
    # _sec_coloc_rpt = '/Users/admin/Downloads/coloc_results/h2004/coloc_output_20221212141510.tsv.gz'
    # _sec_smr_rpt = '/Users/admin/Downloads/coloc_results/h2004/smr_output_20221212142103.tsv.gz'
    # _sec_jlim_rpt = None
    # _sec_fasatenloc_rpt = '/Users/admin/Downloads/coloc_results/h2004/fastenloc_output_20221212153046.tsv.gz'
    # _sec_predixcan_rpt = '/Users/admin/Downloads/coloc_results/h2004/predixcan_output_20221212141515.tsv.gz'
    # _sec_ecaviar_rpt = '/Users/admin/Downloads/coloc_results/h2004/ecaviar_output_20221212141934.tsv.gz'
    # _sec_twas_rpt = '/Users/admin/Downloads/coloc_results/h2004/twas_output_20221223222921.tsv.gz'
    # # H0/H2 causal type: 0:H0; 1:H1; 2:H2 r2<=0.4; 3:H2 0.4<r2<=0.7; 4:H2 0.7<r2<=0.9
    # _sec_causal_type = 2
    # _roc_figure_ensemble_path = '/Users/admin/Downloads/coloc_results/ROC_H1_H2_004.png'
    #
    # ----------------
    # ----------------
    # _sec_coloc_rpt = '/Users/admin/Downloads/coloc_results/h20407/coloc_output_20221212142112.tsv.gz'
    # _sec_smr_rpt = '/Users/admin/Downloads/coloc_results/h20407/smr_output_20221212142748.tsv.gz'
    # _sec_jlim_rpt = None
    # _sec_fasatenloc_rpt = '/Users/admin/Downloads/coloc_results/h20407/fastenloc_output_20221212153049.tsv.gz'
    # _sec_predixcan_rpt = '/Users/admin/Downloads/coloc_results/h20407/predixcan_output_20221212142115.tsv.gz'
    # _sec_ecaviar_rpt = '/Users/admin/Downloads/coloc_results/h20407/ecaviar_output_20221212142707.tsv.gz'
    # _sec_twas_rpt = '/Users/admin/Downloads/coloc_results/h20407/twas_output_20221223223832.tsv.gz'
    # # H0/H2 causal type: 0:H0; 1:H1; 2:H2 r2<=0.4; 3:H2 0.4<r2<=0.7; 4:H2 0.7<r2<=0.9
    # _sec_causal_type = 3
    #
    # _roc_figure_ensemble_path = '/Users/admin/Downloads/coloc_results/ROC_H1_H2_0407.png'
    # ----------------
    # ----------------
    # _sec_coloc_rpt = '/Users/admin/Downloads/coloc_results/h20709/coloc_output_20221212142756.tsv.gz'
    # _sec_smr_rpt = '/Users/admin/Downloads/coloc_results/h20709/smr_output_20221212143433.tsv.gz'
    # _sec_jlim_rpt = None
    # _sec_fasatenloc_rpt = '/Users/admin/Downloads/coloc_results/h20709/fastenloc_output_20221212153052.tsv.gz'
    # _sec_predixcan_rpt = '/Users/admin/Downloads/coloc_results/h20709/predixcan_output_20221212142800.tsv.gz'
    # _sec_ecaviar_rpt = '/Users/admin/Downloads/coloc_results/h20709/ecaviar_output_20221212143404.tsv.gz'
    # _sec_twas_rpt = '/Users/admin/Downloads/coloc_results/h20709/twas_output_20221223224746.tsv.gz'
    # # H0/H2 causal type: 0:H0; 1:H1; 2:H2 r2<=0.4; 3:H2 0.4<r2<=0.7; 4:H2 0.7<r2<=0.9
    # _sec_causal_type = 4
    #
    # _roc_figure_ensemble_path = '/Users/admin/Downloads/coloc_results/ROC_H1_H2_0709.png'
    # ----------------
    #
    # plot_all_against_ensemble_roc(_generated_file,
    #                               _h1_coloc_rpt, _sec_coloc_rpt,
    #                               _h1_smr_rpt, _sec_smr_rpt,
    #                               _h1_jlim_rpt, _sec_jlim_rpt,
    #                               h1_fastenloc_rpt=_h1_fastenloc_rpt, sec_fastenloc_rpt=_sec_fasatenloc_rpt,
    #                               h1_predixcan_rpt=_h1_predixcan_rpt, sec_predixcan_rpt=_sec_predixcan_rpt,
    #                               h1_ecaviar_rpt=_h1_ecaviar_rpt, sec_ecaviar_rpt=_sec_ecaviar_rpt,
    #                               h1_twas_rpt=_h1_twas_rpt, sec_twas_rpt=_sec_twas_rpt,
    #                               sec_causal_type=_sec_causal_type,
    #                               output_figure_path=_roc_figure_ensemble_path, output_dir=_output_dir)
