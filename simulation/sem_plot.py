import os
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

import ranking.intact as intact
import ranking.rra as rra
import sim_result_plot
from ranking.constants import TOOL_SIG_COL_INFO, GENE_ID_COL_NAME, RESULT_TYPE_PVAL, RESULT_TYPE_PROB
from threshold import calc

prob_col_name = 'PROB'
is_positive_col_name = 'IS_POSITIVE'


def calc_stat():
    target_dir = '/Users/admin/Downloads/sim_scheme/split_no_limit'
    ref_file = '/Users/admin/Downloads/sim_scheme/sim_data/gwas_truth.tsv'
    truth_df = pd.read_table(ref_file, usecols=['gene_id', 'causal'])
    truth_df['gene_in_result_tool_cnt'] = 0
    for f in os.listdir(target_dir):
        if not f.startswith('ENSG'):
            continue
        gene_id = f
        gene_dir = os.path.join(target_dir, f)
        tool_cnt = 0
        for f2 in os.listdir(gene_dir):
            f2_full_path = os.path.join(gene_dir, f2)
            if f2.startswith('coloc_output_') or f2.startswith('predixcan_output_') or f2.startswith(
                    'smr_output_') or f2.startswith('fastenloc_output_') or f2.startswith(
                'twas_output_') or f2.startswith('ecaviar_output_'):
                h1_coloc = f2_full_path
                df = pd.read_table(f2_full_path, usecols=['gene_id'])
                df.drop(index=df[~(df['gene_id'] == gene_id)].index, inplace=True)
                if df.shape[0] > 0:
                    tool_cnt += 1
                del df
        truth_df.loc[truth_df['gene_id'] == gene_id, 'gene_in_result_tool_cnt'] = tool_cnt
    truth_df.to_csv(os.path.join(target_dir, 'stat.tsv'), sep='\t', header=True, index=False)


def merge_result():
    target_dir = '/Users/admin/Downloads/sim_scheme/split_no_limit'
    ref_file = '/Users/admin/Downloads/sim_scheme/sim_data/gwas_truth.tsv'
    # truth_df = pd.read_table(ref_file, usecols=['gene_id', 'causal'])
    # truth_df['gene_id'] = truth_df['gene_id'] + '_' + truth_df['gene_id']
    coloc_results = []
    ecaviar_results = []
    smr_results = []
    twas_results = []
    fastenloc_results = []
    predixcan_results = []
    for f in os.listdir(target_dir):
        if not f.startswith('ENSG'):
            continue
        gene_id = f
        gene_dir = os.path.join(target_dir, f)
        for f2 in os.listdir(gene_dir):
            if f2.startswith('.'):
                continue
            f2_full_path = os.path.join(gene_dir, f2)
            df = pd.read_table(f2_full_path)
            # df['gene_id'] = gene_id + '_' + df['gene_id']
            # df.drop(index=df[df['gene_id'] != gene_id].index, inplace=True)
            if df.shape[0] == 0:
                continue
            if f2.startswith('coloc_output_'):
                coloc_results.append(df)
            elif f2.startswith('smr_output_'):
                smr_results.append(df)
            elif f2.startswith('twas_output_'):
                twas_results.append(df)
            elif f2.startswith('ecaviar_output_'):
                ecaviar_results.append(df)
            elif f2.startswith('fastenloc_output_'):
                fastenloc_results.append(df)
            elif f2.startswith('predixcan_output_'):
                predixcan_results.append(df)
    ecaviar_df = pd.concat(ecaviar_results)
    ecaviar_df.sort_values(by='clpp', ascending=False, inplace=True)
    ecaviar_df.drop_duplicates(subset='gene_id', inplace=True)
    ecaviar_df.to_csv(os.path.join(target_dir, 'ecaviar_merged.tsv'), sep='\t', index=False)

    coloc_df = pd.concat(coloc_results)
    coloc_df.sort_values(by=['overall_H4', 'SNP.PP.H4'], ascending=False, inplace=True)
    coloc_df.drop_duplicates(subset='gene_id', inplace=True)
    coloc_df.to_csv(os.path.join(target_dir, 'coloc_merged.tsv'), sep='\t', index=False)

    smr_df = pd.concat(smr_results)
    smr_df.sort_values(by='p_SMR', ascending=True, inplace=True)
    smr_df.drop_duplicates(subset='gene_id', inplace=True)
    smr_df.to_csv(os.path.join(target_dir, 'smr_merged.tsv'), sep='\t', index=False)

    twas_df = pd.concat(twas_results)
    twas_df.sort_values(by='TWAS.P', ascending=True, inplace=True)
    twas_df.drop_duplicates(subset='gene_id', inplace=True)
    twas_df.to_csv(os.path.join(target_dir, 'twas_merged.tsv'), sep='\t', index=False)

    fastenloc_df = pd.concat(fastenloc_results)
    fastenloc_df.sort_values(by='LCP', ascending=False, inplace=True)
    fastenloc_df.drop_duplicates(subset='gene_id', inplace=True)
    fastenloc_df.to_csv(os.path.join(target_dir, 'fastenloc_merged.tsv'), sep='\t', index=False)

    predixcan_df = pd.concat(predixcan_results)
    predixcan_df.sort_values(by='pvalue', ascending=True, inplace=True)
    predixcan_df.drop_duplicates(subset='gene_id', inplace=True)
    predixcan_df.to_csv(os.path.join(target_dir, 'predixcan_merged.tsv'), sep='\t', index=False)


# # refer to sim_result_plot, sem specified method start==========================================
def retrieve_std_df(generated_list):
    std_df = pd.read_table(generated_list, header=None, sep=r'\s+', usecols=[0, 3])
    std_df.columns = [GENE_ID_COL_NAME, is_positive_col_name]
    # std_df[GENE_ID_COL_NAME] = std_df[GENE_ID_COL_NAME] + '_' + std_df[GENE_ID_COL_NAME]
    return std_df


def retrieve_positive_df(generated_list):
    positive_df = pd.read_table(generated_list, header=None, sep=r'\s+', usecols=[0, 3])
    positive_df.columns = [GENE_ID_COL_NAME, is_positive_col_name]
    positive_df.drop(index=positive_df[positive_df[is_positive_col_name] == 0].index, inplace=True)
    positive_df.reset_index(drop=True, inplace=True)
    return positive_df


def prepare_plot_data(generated_list, h1_report,
                      rpt_prob_col_name=None, rpt_pval_col_name=None, tool=None):
    if rpt_prob_col_name is None and rpt_pval_col_name is None:
        raise ValueError('p-value column name and probability column name can not be both null')
    std_df = retrieve_std_df(generated_list)

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
        h1_report_df = pd.read_table(h1_report, usecols=reading_cols)
        h1_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        h1_report_df[prob_col_name] = 1 - h1_report_df[rpt_pval_col_name]
    else:
        h1_report_df = pd.read_table(h1_report, usecols=[rpt_prob_col_name, GENE_ID_COL_NAME])
        h1_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        h1_report_df[prob_col_name] = h1_report_df[rpt_prob_col_name]
    report_df = h1_report_df
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


def plot_all_against_ensemble_prc(generated_file_path, h1_rpt_obj,
                                  output_figure_path=None,
                                  output_dir=''):
    plt.figure().clear()
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    rpt_obj = {}
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            tool_df = prepare_plot_data(generated_file_path, h1_rpt,
                                        rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                        rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                        tool=tool)
            if tool_df.empty:
                continue
            sim_result_plot.plot_prc_curve(tool_df[is_positive_col_name], tool_df[prob_col_name], tool=tool)
            tool_output = os.path.join(output_dir, f'{tool}_result.tsv')
            tool_df.to_csv(tool_output, sep='\t', header=True, index=False)
            rpt_obj[tool] = tool_output
    if len(rpt_obj) == 0:
        print('all report is empty, nothing to do')
        return
    std_df = retrieve_std_df(generated_file_path)
    # --------
    # rra_geo_output_file = os.path.join(output_dir, 'rra_geo.tsv')
    # rra.run_ranking(output_file_path=rra_geo_output_file, rpt=rpt_obj, sample_size=500, method='GEO')
    # rgeo_df = pd.read_table(rra_geo_output_file, usecols=[GENE_ID_COL_NAME, 'geo_p_value'])
    # # Convert rra p_value to probability. TODO this method is not good
    # rgeo_df[prob_col_name] = 1 - rgeo_df['geo_p_value']
    # rgeo_df = pd.merge(left=rgeo_df, right=std_df,
    #                    left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
    #                    how='left')
    # tp_na_bool_series = rgeo_df[is_positive_col_name].isna()
    # rgeo_df.loc[rgeo_df[tp_na_bool_series].index, is_positive_col_name] = 0
    # plot_prc_curve(rgeo_df[is_positive_col_name], rgeo_df[prob_col_name], tool='rGEO')
    # --------intact_expit
    intact_expit_output_file = os.path.join(output_dir, 'intact_expit.tsv')
    intact.run_ranking(rpt=rpt_obj, output_file_path=intact_expit_output_file, prior_fun='expit')
    intact_expit_df = pd.read_table(intact_expit_output_file)
    intact_expit_df = pd.merge(left=intact_expit_df, right=std_df,
                               left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                               how='left')
    tp_na_bool_series = intact_expit_df[is_positive_col_name].isna()
    intact_expit_df.loc[intact_expit_df[tp_na_bool_series].index, is_positive_col_name] = 0
    intact_expit_df.loc[intact_expit_df['intact_probability'].isna(), 'intact_probability'] = 0
    intact_expit_df.to_csv(os.path.join(output_dir, f'intact_expit_merged_std.tsv'), sep='\t', header=True, index=False,
                           na_rep='NA')
    sim_result_plot.plot_prc_curve(intact_expit_df[is_positive_col_name], intact_expit_df['intact_probability'],
                                   tool='intact expit')
    # --------
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


def plot_all_against_ensemble_roc(generated_file_path, h1_rpt_obj,
                                  output_figure_path=None,
                                  output_dir=''):
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
        if h1_rpt is not None and Path(h1_rpt).exists() and os.path.getsize(h1_rpt) > 0:
            tool_df = prepare_plot_data(generated_file_path, h1_rpt,
                                        rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                        rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                        tool=tool)
            if tool_df.empty:
                continue
            sim_result_plot.plot_roc_curve(tool_df[is_positive_col_name], tool_df[prob_col_name], tool=tool)
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
    std_df = retrieve_std_df(generated_file_path)
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
    sim_result_plot.plot_roc_curve(rgeo_df[is_positive_col_name], rgeo_df[prob_col_name], tool='rGEO')
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
    sim_result_plot.plot_roc_curve(intact_expit_df[is_positive_col_name], intact_expit_df['intact_probability'],
                                   tool='intact expit')
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


def plot_precision_recall_f1_comb_bar(
        generated_file_path=None,
        h1_rpt_obj=None,
        output_figure_prefix=None,
        typ='UNION'):
    rpt_obj = {}
    thresh_obj = calc.calc_threshold(h1_rpt_obj)
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        tool_df = prepare_plot_data(generated_file_path, h1_rpt,
                                    rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                    rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                    tool=tool)
        tool_df['mark'] = pd.NA
        if sig_type == RESULT_TYPE_PVAL:
            if tool == 'smr':
                positive_series = (tool_df[sig_column] < thresh_obj[tool]) & (tool_df['p_HEIDI'] > 0.05)
                tool_df['mark'].mask(positive_series, 1, inplace=True)
                tool_df['mark'].mask(~positive_series, 0, inplace=True)
            else:
                tool_df['mark'].mask(tool_df[sig_column] < thresh_obj[tool], 1, inplace=True)
                tool_df['mark'].mask(tool_df[sig_column] >= thresh_obj[tool], 0, inplace=True)
        else:
            tool_df['mark'].mask(tool_df[sig_column] > thresh_obj[tool], 1, inplace=True)
            tool_df['mark'].mask(tool_df[sig_column] <= thresh_obj[tool], 0, inplace=True)
        rpt_obj[tool] = tool_df[[GENE_ID_COL_NAME, 'mark', is_positive_col_name]]

    tools = [tool for tool, _, _ in TOOL_SIG_COL_INFO]
    count = []
    precision_means = []
    precision_stds = []
    recall_means = []
    recall_stds = []
    f1_means = []
    f1_stds = []
    coms = []
    precisions = []
    recalls = []
    f1s = []
    for n in range(1, len(tools) + 1):
        count.append(n)
        ntool_precisions = []
        ntool_recalls = []
        ntool_f1s = []
        n_coms = []
        for com in combinations(tools, n):
            # union n tool results, for the same gene across diff tools,
            # keep the one with the largest mark if typ=UNION,
            # else keep the one with the smallest mark
            df = pd.concat([rpt_obj[tool] for tool in com])
            df.sort_values(by='mark', ascending=typ.upper() != 'UNION', inplace=True)
            df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
            if df.empty:
                ntool_precisions.append(0)
                ntool_recalls.append(0)
                ntool_f1s.append(0)
                continue
            tp = sum((df[f'mark'] == 1) & (df[is_positive_col_name] == 1))
            fp = sum((df[is_positive_col_name] == 0) & (df[f'mark'] == 1))
            fn = sum((df[is_positive_col_name] == 1) & (df[f'mark'] == 0))
            if tp + fp == 0:
                precision = 0
            else:
                precision = tp / (tp + fp)
            ntool_precisions.append(precision)
            if tp + fn == 0:
                recall = 0
            else:
                recall = tp / (tp + fn)
            ntool_recalls.append(recall)
            if 2 * tp + fp + fn == 0:
                f1 = 0
            else:
                f1 = 2 * tp / (2 * tp + fp + fn)
            ntool_f1s.append(f1)
            n_coms.append(com)
        n_union_df = pd.DataFrame({'precision': ntool_precisions, 'recall': ntool_recalls, 'f1': ntool_f1s})
        precision_means.append(n_union_df['precision'].mean())
        precision_stds.append(0 if n_union_df.shape[0] == 1 else n_union_df['precision'].std(ddof=0))
        recall_means.append(n_union_df['recall'].mean())
        recall_stds.append(0 if n_union_df.shape[0] == 1 else n_union_df['recall'].std(ddof=0))
        f1_means.append(n_union_df['f1'].mean())
        f1_stds.append(0 if n_union_df.shape[0] == 1 else n_union_df['f1'].std(ddof=0))
        coms.append(n_coms)
        precisions.append(ntool_precisions)
        recalls.append(ntool_recalls)
        f1s.append(ntool_f1s)

    # merge all genes
    gene_list = []
    for tool in tools:
        gene_list.append(rpt_obj[tool][[GENE_ID_COL_NAME, is_positive_col_name]])
    gene_list_df = pd.concat(gene_list)
    del gene_list
    gene_list_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    # merge all results in rpt_obj
    mark_cols = []
    merged_df = gene_list_df
    for tool in tools:
        merged_df = pd.merge(left=merged_df, right=rpt_obj[tool][[GENE_ID_COL_NAME, 'mark']],
                             left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                             how='outer')
        merged_df.rename(columns={'mark': f'{tool}_mark'}, inplace=True)
        mark_cols.append(f'{tool}_mark')
    merged_df['mvote_positive'] = merged_df[mark_cols].sum(axis=1) >= 3
    mv_tp = merged_df[merged_df['mvote_positive'] & (merged_df[is_positive_col_name] == 1)].shape[0]
    mv_fp = merged_df[(merged_df[is_positive_col_name] == 0) & merged_df['mvote_positive']].shape[0]
    mv_fn = merged_df[(merged_df[is_positive_col_name] == 1) & (~merged_df['mvote_positive'])].shape[0]
    if mv_tp + mv_fp == 0:
        mv_precision = 0
    else:
        mv_precision = mv_tp / (mv_tp + mv_fp)
    if mv_tp + mv_fn == 0:
        mv_recall = 0
    else:
        mv_recall = mv_tp / (mv_tp + mv_fn)
    mv_f1 = 2 * mv_tp / (2 * mv_tp + mv_fp + mv_fn)
    precision_means.append(mv_precision)
    precision_stds.append(0)
    recall_means.append(mv_recall)
    recall_stds.append(0)
    f1_means.append(mv_f1)
    f1_stds.append(0)
    xticklbl = [str(c) for c in count]
    xticklbl.append(f'm_vote\nT=cnt>=3')
    count.append(7)

    # precision mean
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    rects = ax.bar(count, precision_means, 0.5, label='Precision Mean')
    ax.bar_label(rects, fmt='{:0.2f}')
    ax.errorbar(count, precision_means, yerr=precision_stds, fmt=',', ecolor='black', capsize=5)
    ax.set(xticks=[*count], xticklabels=xticklbl)
    ax.set_ylabel('Precision')
    ax.set_xlabel('Num of tools')
    ax.set_title(f'Precision Mean of different combinations of tools')
    plt.savefig(f'{output_figure_prefix}_precision.png')
    plt.close()

    # precision scatter
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    for idx, prs in enumerate(precisions):
        ax.scatter([idx + 1] * len(prs), prs)
        min_precision = min(prs)
        max_precision = max(prs)
        min_precision_com = coms[idx][prs.index(min_precision)]
        max_precision_com = coms[idx][prs.index(max_precision)]
        # for single_prs_idx, single_prs in enumerate(prs):
        #     ax.text(idx + 1, single_prs, ','.join(list(map(lambda t: t[0:1], coms[idx][single_prs_idx]))))
        # show label for min/max only
        ax.text(idx + 1, min_precision, ','.join(list(map(lambda t: t[0:1], min_precision_com))))
        ax.text(idx + 1, max_precision, ','.join(list(map(lambda t: t[0:1], max_precision_com))))
        print(f'for {len(min_precision_com)} tools: '
              f'min_precision: {min_precision}, min_precision_com: {min_precision_com}, '
              f'max_precision: {max_precision}, max_precision_com: {max_precision_com}')
    ax.set_ylabel('Precision')
    ax.set_xlabel('Num of tools')
    ax.set_title(f'Precision Mean of different combinations of tools')
    plt.savefig(f'{output_figure_prefix}_precision_scatter.png')
    plt.close()

    # recall mean
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    rects = ax.bar(count, recall_means, 0.5, label='Recall Mean')
    ax.bar_label(rects, fmt='{:0.2f}')
    ax.errorbar(count, recall_means, yerr=recall_stds, fmt=',', ecolor='black', capsize=5)
    ax.set(xticks=[*count], xticklabels=xticklbl)
    ax.set_ylabel('Recall')
    ax.set_xlabel('Num of tools')
    ax.set_title(f'Recall Mean of different combinations of tools')
    plt.savefig(f'{output_figure_prefix}_recall.png')
    plt.close()

    # recall scatter
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    for idx, rcl in enumerate(recalls):
        ax.scatter([idx + 1] * len(rcl), rcl)
        min_recall = min(rcl)
        max_recall = max(rcl)
        min_recall_com = coms[idx][rcl.index(min_recall)]
        max_recall_com = coms[idx][rcl.index(max_recall)]
        # for single_rcl_idx, single_rcl in enumerate(rcl):
        #     ax.text(idx + 1, single_rcl, ','.join(list(map(lambda t: t[0:1], coms[idx][single_rcl_idx]))))
        # show label for min/max only
        ax.text(idx + 1, min_recall, ','.join(list(map(lambda t: t[0:1], min_recall_com))))
        ax.text(idx + 1, max_recall, ','.join(list(map(lambda t: t[0:1], max_recall_com))))
        print(f'for {len(min_recall_com)} tools: '
              f'min_recall: {min_recall}, min_recall_com: {min_recall_com}, '
              f'max_recall: {max_recall}, max_recall_com: {max_recall_com}')
    ax.set_ylabel('recall')
    ax.set_xlabel('Num of tools')
    ax.set_title(f'recall Mean of different combinations of tools')
    plt.savefig(f'{output_figure_prefix}_recall_scatter.png')
    plt.close()

    # f1 mean
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    rects = ax.bar(count, f1_means, 0.5, label='F1 Mean')
    ax.bar_label(rects, fmt='{:0.2f}')
    ax.errorbar(count, f1_means, yerr=f1_stds, fmt=',', ecolor='black', capsize=5)
    ax.set(xticks=[*count], xticklabels=xticklbl)
    ax.set_ylabel('F1')
    ax.set_xlabel('Num of tools')
    ax.set_title(f'F1 Mean of different combinations of tools')
    plt.savefig(f'{output_figure_prefix}_f1.png')
    plt.close()

    # f1 scatter
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    for idx, cf1 in enumerate(f1s):
        ax.scatter([idx + 1] * len(cf1), cf1)
        min_f1 = min(cf1)
        max_f1 = max(cf1)
        min_f1_com = coms[idx][cf1.index(min_f1)]
        max_f1_com = coms[idx][cf1.index(max_f1)]
        # for single_f1_idx, single_f1 in enumerate(rcl):
        #     ax.text(idx + 1, single_f1, ','.join(list(map(lambda t: t[0:1], coms[idx][single_f1_idx]))))
        # show label for min/max only
        ax.text(idx + 1, min_f1, ','.join(list(map(lambda t: t[0:1], min_f1_com))))
        ax.text(idx + 1, max_f1, ','.join(list(map(lambda t: t[0:1], max_f1_com))))
        print(f'for {len(min_f1_com)} tools: '
              f'min_f1: {min_f1}, min_f1_com: {min_f1_com}, '
              f'max_f1: {max_f1}, max_f1_com: {max_f1_com}')
    ax.set_ylabel('F1')
    ax.set_xlabel('Num of tools')
    ax.set_title(f'F1 Mean of different combinations of tools')
    plt.savefig(f'{output_figure_prefix}_F1_scatter.png')
    plt.close()


# # refer to sim_result_plot, sem specified method end==========================================

def plot_tools_compare_precision_recall_f1_bar(
        generated_file_path=None,
        h1_rpt_obj=None,
        output_figure_prefix=None):
    rpt_obj = {}
    thresh_obj = calc.calc_threshold(h1_rpt_obj)
    tools = [tool for tool, _, _ in TOOL_SIG_COL_INFO]
    xticklbl = [tool for tool in tools]
    precisions = []
    recalls = []
    f1s = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        h1_rpt = h1_rpt_obj.get(tool)
        tool_df = prepare_plot_data(generated_file_path, h1_rpt,
                                    rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                    rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                    tool=tool)
        tool_df['mark'] = pd.NA
        if sig_type == RESULT_TYPE_PVAL:
            if tool == 'smr':
                positive_series = (tool_df[sig_column] < thresh_obj[tool]) & (tool_df['p_HEIDI'] > 0.05)
                tool_df['mark'].mask(positive_series, 1, inplace=True)
                tool_df['mark'].mask(~positive_series, 0, inplace=True)
            else:
                tool_df['mark'].mask(tool_df[sig_column] < thresh_obj[tool], 1, inplace=True)
                tool_df['mark'].mask(tool_df[sig_column] >= thresh_obj[tool], 0, inplace=True)
        else:
            tool_df['mark'].mask(tool_df[sig_column] > thresh_obj[tool], 1, inplace=True)
            tool_df['mark'].mask(tool_df[sig_column] <= thresh_obj[tool], 0, inplace=True)
        rpt_obj[tool] = tool_df[[GENE_ID_COL_NAME, 'mark', is_positive_col_name]]
        # calc tool precision/recall/f1
        tp = sum((tool_df[f'mark'] == 1) & (tool_df[is_positive_col_name] == 1))
        fp = sum((tool_df[is_positive_col_name] == 0) & (tool_df[f'mark'] == 1))
        fn = sum((tool_df[is_positive_col_name] == 1) & (tool_df[f'mark'] == 0))
        if tp + fp == 0:
            precision = 0
        else:
            precision = tp / (tp + fp)
        precisions.append(precision)
        if tp + fn == 0:
            recall = 0
        else:
            recall = tp / (tp + fn)
        recalls.append(recall)
        if 2 * tp + fp + fn == 0:
            f1 = 0
        else:
            f1 = 2 * tp / (2 * tp + fp + fn)
        f1s.append(f1)
    # calc majority vote precision/recall/f1
    gene_list = []
    for tool in tools:
        gene_list.append(rpt_obj[tool][[GENE_ID_COL_NAME, is_positive_col_name]])
    gene_list_df = pd.concat(gene_list)
    del gene_list
    gene_list_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    # merge all results in rpt_obj
    mark_cols = []
    merged_df = gene_list_df
    for tool in tools:
        merged_df = pd.merge(left=merged_df, right=rpt_obj[tool][[GENE_ID_COL_NAME, 'mark']],
                             left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                             how='outer')
        merged_df.rename(columns={'mark': f'{tool}_mark'}, inplace=True)
        mark_cols.append(f'{tool}_mark')
    merged_df['mvote_positive'] = merged_df[mark_cols].sum(axis=1) >= 3
    mv_tp = merged_df[merged_df['mvote_positive'] & (merged_df[is_positive_col_name] == 1)].shape[0]
    mv_fp = merged_df[(merged_df[is_positive_col_name] == 0) & merged_df['mvote_positive']].shape[0]
    mv_fn = merged_df[(merged_df[is_positive_col_name] == 1) & (~merged_df['mvote_positive'])].shape[0]
    if mv_tp + mv_fp == 0:
        mv_precision = 0
    else:
        mv_precision = mv_tp / (mv_tp + mv_fp)
    if mv_tp + mv_fn == 0:
        mv_recall = 0
    else:
        mv_recall = mv_tp / (mv_tp + mv_fn)
    mv_f1 = 2 * mv_tp / (2 * mv_tp + mv_fp + mv_fn)

    precisions.append(mv_precision)
    recalls.append(mv_recall)
    f1s.append(mv_f1)
    xticklbl.append(f'm_vote\nT=cnt>=3')

    # precision
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    rects = ax.bar([*range(1, len(xticklbl) + 1)], precisions, 0.5, label='Precision')
    ax.bar_label(rects, fmt='{:0.2f}')
    ax.set(xticks=[*range(1, len(xticklbl) + 1)], xticklabels=xticklbl)
    ax.set_ylabel('Precision')
    ax.set_xlabel('Tools')
    ax.set_title(f'Precision different tools')
    plt.savefig(f'{output_figure_prefix}_precision.png')
    plt.close()

    # recall
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    rects = ax.bar([*range(1, len(xticklbl) + 1)], recalls, 0.5, label='Recall')
    ax.bar_label(rects, fmt='{:0.2f}')
    ax.set(xticks=[*range(1, len(xticklbl) + 1)], xticklabels=xticklbl)
    ax.set_ylabel('Recall')
    ax.set_xlabel('Tools')
    ax.set_title(f'Recall of different tools')
    plt.savefig(f'{output_figure_prefix}_recall.png')
    plt.close()

    # f1
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(8, 6))
    rects = ax.bar([*range(1, len(xticklbl) + 1)], f1s, 0.5, label='F1')
    ax.bar_label(rects, fmt='{:0.2f}')
    ax.set(xticks=[*range(1, len(xticklbl) + 1)], xticklabels=xticklbl)
    ax.set_ylabel('F1')
    ax.set_xlabel('Tools')
    ax.set_title(f'F1 of different tools')
    plt.savefig(f'{output_figure_prefix}_f1.png')
    plt.close()


def download_from_28(file_list_file, dest_dir):
    with open(file_list_file) as f_list:
        for f in f_list:
            f = f.strip()
            splits = f.split(os.sep)
            trait = splits[-6]
            file_name = splits[-1]
            trait_dir = os.path.join(dest_dir, trait)
            Path(trait_dir).mkdir(exist_ok=True, parents=True)
            dest_file = os.path.join(trait_dir, file_name)
            print(f'Downloading {f} to dest {dest_file}')
            os.system(f'cp {f} {dest_file}')
            if not os.path.exists(dest_file):
                raise ValueError(f'{f} synced failed')


def plot():
    # # ==========================================
    # prc_path = os.path.join('/Users/admin/Downloads/sim_scheme/split_no_limit',
    #                         f'PRC.png')
    # plot_all_against_ensemble_prc(
    #     '/Users/admin/Downloads/sim_scheme/sim_data/gwas.truth',
    #     {
    #         'twas': '/Users/admin/Downloads/sim_scheme/split_no_limit/twas_merged.tsv',
    #         'smr': '/Users/admin/Downloads/sim_scheme/split_no_limit/smr_merged.tsv',
    #         'coloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/coloc_merged.tsv',
    #         'ecaviar': '/Users/admin/Downloads/sim_scheme/split_no_limit/ecaviar_merged.tsv',
    #         'fastenloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/fastenloc_merged.tsv',
    #         'predixcan': '/Users/admin/Downloads/sim_scheme/split_no_limit/predixcan_merged.tsv'},
    #     output_figure_path=prc_path,
    #     output_dir='/Users/admin/Downloads/sim_scheme/split_no_limit/')
    #
    # roc_path = os.path.join('/Users/admin/Downloads/sim_scheme/split_no_limit',
    #                         f'ROC.png')
    # plot_all_against_ensemble_roc(
    #     '/Users/admin/Downloads/sim_scheme/sim_data/gwas.truth',
    #     {
    #         'twas': '/Users/admin/Downloads/sim_scheme/split_no_limit/twas_merged.tsv',
    #         'smr': '/Users/admin/Downloads/sim_scheme/split_no_limit/smr_merged.tsv',
    #         'coloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/coloc_merged.tsv',
    #         'ecaviar': '/Users/admin/Downloads/sim_scheme/split_no_limit/ecaviar_merged.tsv',
    #         'fastenloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/fastenloc_merged.tsv',
    #         'predixcan': '/Users/admin/Downloads/sim_scheme/split_no_limit/predixcan_merged.tsv'},
    #     output_figure_path=roc_path,
    #     output_dir='/Users/admin/Downloads/sim_scheme/split_no_limit/')
    # ==========================================
    # prefix_path = os.path.join('/Users/admin/Downloads/sim_scheme/split_no_limit',
    #                            f'UNION')
    # plot_precision_recall_f1_comb_bar(
    #     '/Users/admin/Downloads/sim_scheme/sim_data/gwas.truth',
    #     {
    #         'twas': '/Users/admin/Downloads/sim_scheme/split_no_limit/twas_merged.tsv',
    #         'smr': '/Users/admin/Downloads/sim_scheme/split_no_limit/smr_merged.tsv',
    #         'coloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/coloc_merged.tsv',
    #         'ecaviar': '/Users/admin/Downloads/sim_scheme/split_no_limit/ecaviar_merged.tsv',
    #         'fastenloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/fastenloc_merged.tsv',
    #         'predixcan': '/Users/admin/Downloads/sim_scheme/split_no_limit/predixcan_merged.tsv'
    #     },
    #     output_figure_prefix=prefix_path,
    #     typ='UNION')
    #
    # prefix_path2 = os.path.join('/Users/admin/Downloads/sim_scheme/split_no_limit',
    #                             f'INTER')
    # plot_precision_recall_f1_comb_bar(
    #     '/Users/admin/Downloads/sim_scheme/sim_data/gwas.truth',
    #     {
    #         'twas': '/Users/admin/Downloads/sim_scheme/split_no_limit/twas_merged.tsv',
    #         'smr': '/Users/admin/Downloads/sim_scheme/split_no_limit/smr_merged.tsv',
    #         'coloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/coloc_merged.tsv',
    #         'ecaviar': '/Users/admin/Downloads/sim_scheme/split_no_limit/ecaviar_merged.tsv',
    #         'fastenloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/fastenloc_merged.tsv',
    #         'predixcan': '/Users/admin/Downloads/sim_scheme/split_no_limit/predixcan_merged.tsv'
    #     },
    #     output_figure_prefix=prefix_path2,
    #     typ='INTER')
    #
    # heat_path = '/Users/admin/Downloads/sim_scheme/split_no_limit/KeepNA.png'
    # sim_result_plot.plot_spearman_heatmap({
    #     'twas': '/Users/admin/Downloads/sim_scheme/split_no_limit/twas_merged.tsv',
    #     'smr': '/Users/admin/Downloads/sim_scheme/original_smr_no_p_limit/smr_merged.tsv',
    #     'coloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/coloc_merged.tsv',
    #     'ecaviar': '/Users/admin/Downloads/sim_scheme/split_no_limit/ecaviar_merged.tsv',
    #     'fastenloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/fastenloc_merged.tsv',
    #     'predixcan': '/Users/admin/Downloads/sim_scheme/split_no_limit/predixcan_merged.tsv'
    # }, tested_gene_df=retrieve_std_df('/Users/admin/Downloads/sim_scheme/sim_data/gwas.truth'),
    #     output_figure_path=heat_path)

    prefix_compare_plot = os.path.join('/Users/admin/Downloads/sim_scheme/split_no_limit',
                                       'compare')
    plot_tools_compare_precision_recall_f1_bar(
        '/Users/admin/Downloads/sim_scheme/sim_data/gwas.truth',
        {
            'twas': '/Users/admin/Downloads/sim_scheme/split_no_limit/twas_merged.tsv',
            'smr': '/Users/admin/Downloads/sim_scheme/split_no_limit/smr_merged.tsv',
            'coloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/coloc_merged.tsv',
            'ecaviar': '/Users/admin/Downloads/sim_scheme/split_no_limit/ecaviar_merged.tsv',
            'fastenloc': '/Users/admin/Downloads/sim_scheme/split_no_limit/fastenloc_merged.tsv',
            'predixcan': '/Users/admin/Downloads/sim_scheme/split_no_limit/predixcan_merged.tsv'
        },
        output_figure_prefix=prefix_compare_plot)


if __name__ == '__main__':
    # download_from_28('/Users/admin/Downloads/sim_scheme/split_no_limit/files.txt',
    #                  '/Users/admin/Downloads/sim_scheme/split_no_limit')
    # merge_result()
    plot()
