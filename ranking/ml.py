import os
from pathlib import Path
import sys
import joblib
import pandas as pd
import xgboost as xgb
from sklearn import ensemble
from sklearn import linear_model
from sklearn import metrics
from sklearn import svm
from sklearn.model_selection import train_test_split
from pulearn import BaggingPuClassifier

from ranking.constants import GENE_ID_COL_NAME
from ranking.constants import TOOL_SIG_COL_INFO
from ranking.constants import RESULT_TYPE_PVAL
from ranking.constants import RESULT_TYPE_PROB
from ranking.constants import AVG_RANKING_COL_NAME

PREDICT_COL_NAME = 'ml_probability'
LABEL_COL_NAME = 'label'
SIM_MODEL_ORDER = ['h0', 'h1', 'h2004', 'h20407', 'h20709']


def retrieve_sim_label_df(generated_list, causal_types=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if causal_types is None:
        causal_types = [i for i in range(0, len(SIM_MODEL_ORDER))]
    generated_list_df = pd.read_table(generated_list,
                                      usecols=['gene', 'gwas_causal_snp', 'gwas_max_assoc_p', 'gwas_r2'] +
                                              [f'eqtl_max_assoc_p_{causal_type}' for causal_type in causal_types])
    df_list = []
    for idx, col in enumerate([f'eqtl_max_assoc_p_{causal_type}' for causal_type in causal_types]):
        df = generated_list_df[(generated_list_df['gwas_max_assoc_p'] <= 1.0E-5) & (
                generated_list_df['gwas_r2'] >= 0.8) & (generated_list_df[col] <= 0.01)].copy()
        # h1 label is 1, rest are 0
        df[LABEL_COL_NAME] = 1 if idx == 1 else 0
        df[GENE_ID_COL_NAME] = df['gene'] + f'_{idx}'
        df_list.append(df[[GENE_ID_COL_NAME, LABEL_COL_NAME]])
    std_df = pd.concat(df_list)
    # std_df has 2 columns: [GENE_ID_COL_NAME, LABEL_COL_NAME]
    return std_df


def read_sim_result(tool, rpts, sig_column):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    tool_df_list = []
    for idx, sim_model in enumerate(SIM_MODEL_ORDER):
        report = rpts.get(sim_model)
        if report is None or len(report) == 0 or (not os.path.exists(report)) or os.path.getsize(report) <= 0:
            continue
        report_df = pd.read_table(report, usecols=[sig_column, GENE_ID_COL_NAME])
        report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        report_df[GENE_ID_COL_NAME] = report_df[GENE_ID_COL_NAME] + f'_{idx}'
        tool_df_list.append(report_df)
    result_df = pd.concat(tool_df_list)
    result_df.rename(columns={sig_column: tool}, inplace=True)
    # result_df has 2 columns: [GENE_ID_COL_NAME, tool], tool col means tool result value
    return result_df


def prepare_sim_train_data(generated_file_path, rpts):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    tool_df_list = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        tool_rpts = rpts.get(tool)
        if tool_rpts is None or len(tool_rpts) == 0:
            continue
        tool_df_list.append(read_sim_result(tool, tool_rpts, sig_column))
    ranking_df = None
    for tool_df in tool_df_list:
        if tool_df is None or tool_df.shape[0] == 0:
            continue
        if ranking_df is None:
            ranking_df = tool_df
        else:
            ranking_df = pd.merge(left=ranking_df, right=tool_df,
                                  left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                                  how='outer')
    std_df = retrieve_sim_label_df(generated_file_path)
    result_df = pd.merge(left=std_df, right=ranking_df,
                         left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                         how='outer')
    # label na 置为0
    tp_na_bool_series = result_df[LABEL_COL_NAME].isna()
    result_df.loc[result_df[tp_na_bool_series].index, LABEL_COL_NAME] = 0
    # NOTE: 如果找到的基因集合有可能超过std_df中基因集合(比如silver std数据?), 所以补充缺失的工具数据的操作需要再一个新的for循环中操作
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        if tool not in result_df.columns:
            # NOTE: 预测输入的工具个数和训练时输入工具的个数必须保持一致, 否则报错. 所以需要补充没有结果的工具的数据.
            # 补充缺失工具的数据
            result_df[tool] = 0 if sig_type == RESULT_TYPE_PROB else 1
        else:
            # 有结果的工具补充缺失基因的数据
            # probability na 置为0, p-value na置为1
            if sig_type == RESULT_TYPE_PVAL:
                sig_na_bool_series = result_df[tool].isna()
                result_df.loc[result_df[sig_na_bool_series].index, tool] = 1
            else:
                sig_na_bool_series = result_df[tool].isna()
                result_df.loc[result_df[sig_na_bool_series].index, tool] = 0
    result_df = result_df.reindex(
        columns=[GENE_ID_COL_NAME, LABEL_COL_NAME] + [tool for tool, _, _ in TOOL_SIG_COL_INFO], copy=False)
    return result_df


def train_xgb(x_train, y_train, model_save_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    model = xgb.XGBRegressor()
    model.fit(x_train, y_train, verbose=True)
    if model_save_path is not None:
        model.save_model(model_save_path)
    return model


def train_model_xgb(input_df, test_result_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if input_df is None or input_df.shape[0] == 0:
        print('input is empty, nothing to do')
        return
    x_df = input_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    y_df = input_df[LABEL_COL_NAME]

    x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.2)
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble.ubj')
    if Path(model_path).exists():
        model = xgb.XGBRegressor()
        model.load_model(model_path)
    else:
        model = train_xgb(x_train, y_train, model_save_path=model_path)
    merged_df = pd.DataFrame(y_test)
    merged_df.reset_index(drop=True, inplace=True)
    y_pred = model.predict(x_test)
    merged_df['prediction'] = y_pred
    predictions = [round(value) for value in y_pred]
    print("Accuracy of xgboost:", metrics.accuracy_score(y_test, predictions))
    print("MAE of xgboost:", metrics.mean_absolute_error(y_test, y_pred))
    print("MSE of xgboost:", metrics.mean_squared_error(y_test, y_pred))
    print("R2 of xgboost:", metrics.r2_score(y_test, y_pred))
    merged_df.to_csv(test_result_file, sep='\t', header=True, index=False)


def train_xgb_rf(x_train, y_train, model_save_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    model = xgb.XGBRFRegressor()
    model.fit(x_train, y_train, verbose=True)
    if model_save_path is not None:
        model.save_model(model_save_path)
    return model


def train_model_xgb_rf(input_df, test_result_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if input_df is None or input_df.shape[0] == 0:
        print('input is empty, nothing to do')
        return
    x_df = input_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    y_df = input_df[LABEL_COL_NAME]

    x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.2)
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_rf.ubj')
    if Path(model_path).exists():
        model = xgb.XGBRFRegressor()
        model.load_model(model_path)
    else:
        model = train_xgb_rf(x_train, y_train, model_save_path=model_path)
    merged_df = pd.DataFrame(y_test)
    merged_df.reset_index(drop=True, inplace=True)
    y_pred = model.predict(x_test)
    merged_df['prediction'] = y_pred
    predictions = [round(value) for value in y_pred]
    print("Accuracy of xgboost rf:", metrics.accuracy_score(y_test, predictions))
    print("MAE of xgboost random forest:", metrics.mean_absolute_error(y_test, y_pred))
    print("MSE of xgboost random forest:", metrics.mean_squared_error(y_test, y_pred))
    print("R2 of xgboost random forest:", metrics.r2_score(y_test, y_pred))
    merged_df.to_csv(test_result_file, sep='\t', header=True, index=False)


def train_svm(x_train, y_train, model_save_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    model = svm.LinearSVR(max_iter=10000)
    model.fit(x_train, y_train)
    if model_save_path is not None:
        joblib.dump(model, model_save_path)
    return model


def train_svm_model(input_df, test_result_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if input_df is None or input_df.shape[0] == 0:
        print('input is empty, nothing to do')
        return
    x_df = input_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    y_df = input_df[LABEL_COL_NAME]

    x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.2)
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_svm')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        model = train_svm(x_train, y_train, model_save_path=model_path)
    merged_df = pd.DataFrame(y_test)
    merged_df.reset_index(drop=True, inplace=True)
    y_pred = model.predict(x_test)
    merged_df['prediction'] = y_pred
    predictions = [round(value) for value in y_pred]
    print("Accuracy of svm:", metrics.accuracy_score(y_test, predictions))
    print("MAE of svm:", metrics.mean_absolute_error(y_test, y_pred))
    print("MSE of svm:", metrics.mean_squared_error(y_test, y_pred))
    print("R2 of svm:", metrics.r2_score(y_test, y_pred))
    merged_df.to_csv(test_result_file, sep='\t', header=True, index=False)


def train_rf(x_train, y_train, model_save_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    model = ensemble.RandomForestRegressor(n_estimators=10000, max_depth=10)
    model.fit(x_train, y_train)
    if model_save_path is not None:
        joblib.dump(model, model_save_path)
    return model


def train_rf_model(input_df, test_result_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if input_df is None or input_df.shape[0] == 0:
        print('input is empty, nothing to do')
        return
    x_df = input_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    y_df = input_df[LABEL_COL_NAME]

    x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.2)
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_rf')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        model = train_rf(x_train, y_train, model_save_path=model_path)
    merged_df = pd.DataFrame(y_test)
    merged_df.reset_index(drop=True, inplace=True)
    y_pred = model.predict(x_test)
    merged_df['prediction'] = y_pred
    predictions = [round(value) for value in y_pred]
    print("Accuracy of random forest:", metrics.accuracy_score(y_test, predictions))
    print("MAE of random forest:", metrics.mean_absolute_error(y_test, y_pred))
    print("MSE of random forest:", metrics.mean_squared_error(y_test, y_pred))
    print("R2 of random forest:", metrics.r2_score(y_test, y_pred))
    if test_result_file is not None:
        merged_df.to_csv(test_result_file, sep='\t', header=True, index=False)


def train_lr(x_train, y_train, model_save_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    model = linear_model.LogisticRegression(max_iter=10000)
    model.fit(x_train, y_train)
    if model_save_path is not None:
        joblib.dump(model, model_save_path)
    return model


def train_lr_model(input_df, test_result_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if input_df is None or input_df.shape[0] == 0:
        print('input is empty, nothing to do')
        return
    x_df = input_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    y_df = input_df[LABEL_COL_NAME]

    x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.2)
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_lr')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        model = train_lr(x_train, y_train, model_save_path=model_path)
    merged_df = pd.DataFrame(y_test)
    merged_df.reset_index(drop=True, inplace=True)
    y_pred = model.predict(x_test)
    merged_df['prediction'] = y_pred
    predictions = [round(value) for value in y_pred]
    print("Accuracy of logistic regression:", metrics.accuracy_score(y_test, predictions))
    print("MAE of logistic regression:", metrics.mean_absolute_error(y_test, y_pred))
    print("MSE of logistic regression:", metrics.mean_squared_error(y_test, y_pred))
    print("R2 of logistic regression:", metrics.r2_score(y_test, y_pred))
    merged_df.to_csv(test_result_file, sep='\t', header=True, index=False)


def train_pu(x_train, y_train, model_save_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    rf = ensemble.RandomForestClassifier(n_estimators=10000, max_depth=10)
    pu_estimator = BaggingPuClassifier(
        base_estimator=rf, n_estimators=10000)
    model = pu_estimator.fit(x_train, y_train)
    if model_save_path is not None:
        joblib.dump(model, model_save_path)
    return model


def train_pu_model(input_df, test_result_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if input_df is None or input_df.shape[0] == 0:
        print('input is empty, nothing to do')
        return
    x_df = input_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    y_df = input_df[LABEL_COL_NAME]

    x_train, x_test, y_train, y_test = train_test_split(x_df, y_df, test_size=0.2)
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_pu')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        model = train_pu(x_train, y_train, model_save_path=model_path)
    merged_df = pd.DataFrame(y_test)
    merged_df.reset_index(drop=True, inplace=True)
    y_pred = model.predict(x_test)
    merged_df['prediction'] = y_pred
    predictions = [round(value) for value in y_pred]
    print("Accuracy of PU learning:", metrics.accuracy_score(y_test, predictions))
    print("MAE of PU learning:", metrics.mean_absolute_error(y_test, y_pred))
    print("MSE of PU learning:", metrics.mean_squared_error(y_test, y_pred))
    print("R2 of PU learning:", metrics.r2_score(y_test, y_pred))
    merged_df.to_csv(test_result_file, sep='\t', header=True, index=False)


def read_tool_result(rpt, tool_name, sig_col_name):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None or (not os.path.exists(rpt)) or os.path.getsize(rpt) <= 0:
        return None
    rpt_df = pd.read_table(rpt, usecols=[GENE_ID_COL_NAME, sig_col_name])
    rpt_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    rpt_df.rename(columns={sig_col_name: tool_name}, inplace=True)
    return rpt_df


def prepare_ranking_input(rpts):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpts is None:
        print('No reports provided')
        return None
    non_empty_rpt_cnt = 0
    for rpt in rpts.values():
        non_empty_rpt_cnt += 0 if rpt is None else 1
    if non_empty_rpt_cnt < 2:
        print('At least 2 reports required to run ensemble ranking, now only got {non_empty_rpt_cnt}, nothing to do')
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
    tool_ranking_cols = []
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        if tool not in ranking_df.columns:
            # NOTE: 预测输入的工具个数和训练时输入工具的个数必须保持一致, 否则报错. 所以需要补充没有结果的工具的数据.
            # 补充缺失工具的数据
            ranking_df[tool] = 0 if sig_type == RESULT_TYPE_PROB else 1
        else:
            tool_ranking_cols.append(f'{tool}_ranking')
            # 有结果的工具补充缺失基因的数据
            # probability na 置为0, p-value na置为1
            ranking_df.sort_values(tool, ascending=sig_type == RESULT_TYPE_PVAL, inplace=True)
            sig_na_bool_series = ranking_df[tool].isna()
            ranking_df.loc[ranking_df[sig_na_bool_series].index, tool] = 1 if sig_type == RESULT_TYPE_PVAL else 0
            ranking_df.loc[ranking_df[~sig_na_bool_series].index, f'{tool}_ranking'] = \
                range(1, ranking_df[~sig_na_bool_series].shape[0] + 1)
            ranking_df.loc[ranking_df[sig_na_bool_series].index, f'{tool}_ranking'] = ranking_df.shape[0] + 1

    ranking_df[AVG_RANKING_COL_NAME] = ranking_df[tool_ranking_cols].sum(axis=1) / len(tool_ranking_cols)
    ranking_df.drop(columns=tool_ranking_cols, inplace=True)
    # tool order is important
    prediction_required_cols = [GENE_ID_COL_NAME] + [tool for tool, _, _ in TOOL_SIG_COL_INFO]
    ranking_df = ranking_df.reindex(
        columns=prediction_required_cols + [col for col in ranking_df.columns if col not in prediction_required_cols],
        copy=False)
    return ranking_df


def run_ranking(rpt=None, output_file_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None:
        print('No data provided')
        return None
    target_df = prepare_ranking_input(rpt)
    if target_df is None or target_df.shape[0] == 0:
        print('Data provided is empty')
        return None
    target_df.reset_index(drop=True, inplace=True)
    # missing tool columns will be filled in prepare_ranking_input
    x_df = target_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble.ubj')
    if Path(model_path).exists():
        model = xgb.XGBRegressor()
        model.load_model(model_path)
    else:
        raise ValueError('Model does not exist!')
    target_df[PREDICT_COL_NAME] = model.predict(x_df)
    target_df.sort_values(PREDICT_COL_NAME, ascending=False, inplace=True)
    target_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


def run_ranking_xgb_rf(rpt=None, output_file_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None:
        print('No data provided')
        return None
    target_df = prepare_ranking_input(rpt)
    if target_df is None or target_df.shape[0] == 0:
        print('Data provided is empty')
        return None
    target_df.reset_index(drop=True, inplace=True)
    # missing tool columns will be filled in prepare_ranking_input
    x_df = target_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_rf.ubj')
    if Path(model_path).exists():
        model = xgb.XGBRFRegressor()
        model.load_model(model_path)
    else:
        raise ValueError('Model does not exist!')
    target_df[PREDICT_COL_NAME] = model.predict(x_df)
    target_df.sort_values(PREDICT_COL_NAME, ascending=False, inplace=True)
    target_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


def run_ranking_svm(rpt=None, output_file_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None:
        print('No data provided')
        return None
    target_df = prepare_ranking_input(rpt)
    if target_df is None or target_df.shape[0] == 0:
        print('Data provided is empty')
        return None
    target_df.reset_index(drop=True, inplace=True)
    # missing tool columns will be filled in prepare_ranking_input
    x_df = target_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_svm')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        raise ValueError('Model does not exist!')
    target_df[PREDICT_COL_NAME] = model.predict(x_df)
    target_df.sort_values(PREDICT_COL_NAME, ascending=False, inplace=True)
    target_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


def run_ranking_rf(rpt=None, output_file_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None:
        print('No data provided')
        return None
    target_df = prepare_ranking_input(rpt)
    if target_df is None or target_df.shape[0] == 0:
        print('Data provided is empty')
        return None
    target_df.reset_index(drop=True, inplace=True)
    # missing tool columns will be filled in prepare_ranking_input
    x_df = target_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_rf')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        raise ValueError('Model does not exist!')
    target_df[PREDICT_COL_NAME] = model.predict(x_df)
    target_df.sort_values(PREDICT_COL_NAME, ascending=False, inplace=True)
    target_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


def run_ranking_lr(rpt=None, output_file_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None:
        print('No data provided')
        return None
    target_df = prepare_ranking_input(rpt)
    if target_df is None or target_df.shape[0] == 0:
        print('Data provided is empty')
        return None
    target_df.reset_index(drop=True, inplace=True)
    # missing tool columns will be filled in prepare_ranking_input
    x_df = target_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_lr')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        raise ValueError('Model does not exist!')
    target_df[PREDICT_COL_NAME] = model.predict(x_df)
    target_df.sort_values(PREDICT_COL_NAME, ascending=False, inplace=True)
    target_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


def run_ranking_pu(rpt=None, output_file_path=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if rpt is None:
        print('No data provided')
        return None
    target_df = prepare_ranking_input(rpt)
    if target_df is None or target_df.shape[0] == 0:
        print('Data provided is empty')
        return None
    target_df.reset_index(drop=True, inplace=True)
    # missing tool columns will be filled in prepare_ranking_input
    x_df = target_df[[tool for tool, _, _ in TOOL_SIG_COL_INFO]]
    model_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'model', 'ensemble_pu')
    if Path(model_path).exists():
        model = joblib.load(model_path)
    else:
        raise ValueError('Model does not exist!')
    target_df[PREDICT_COL_NAME] = model.predict(x_df)
    target_df.sort_values(PREDICT_COL_NAME, ascending=False, inplace=True)
    target_df.to_csv(output_file_path, sep='\t', header=True, index=False)
    return output_file_path


if __name__ == '__main__':
    sim_file_set = [
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.1_hsq0.05/output_20000_20000_20221102021042/generated_20221102021042.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.05/sim_tissue_20221102021042_h0/simu_gwas_20221102021042',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.05/sim_tissue_20221102021042_h1/simu_gwas_20221102021042',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.05/sim_tissue_20221102021042_h2004/simu_gwas_20221102021042',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.05/sim_tissue_20221102021042_h20407/simu_gwas_20221102021042',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.05/sim_tissue_20221102021042_h20709/simu_gwas_20221102021042',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.1_hsq0.1/output_20000_20000_20221101170205/generated_20221101170205.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.1/sim_tissue_20221101170205_h0/simu_gwas_20221101170205',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.1/sim_tissue_20221101170205_h1/simu_gwas_20221101170205',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.1/sim_tissue_20221101170205_h2004/simu_gwas_20221101170205',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.1/sim_tissue_20221101170205_h20407/simu_gwas_20221101170205',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.1/sim_tissue_20221101170205_h20709/simu_gwas_20221101170205',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.1_hsq0.2/output_20000_20000_20221102031352/generated_20221102031352.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.2/sim_tissue_20221102031352_h0/simu_gwas_20221102031352',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.2/sim_tissue_20221102031352_h1/simu_gwas_20221102031352',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.2/sim_tissue_20221102031352_h2004/simu_gwas_20221102031352',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.2/sim_tissue_20221102031352_h20407/simu_gwas_20221102031352',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.1_hsq0.2/sim_tissue_20221102031352_h20709/simu_gwas_20221102031352',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.2_hsq0.05/output_20000_20000_20221101115128/generated_20221101115128.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.05/sim_tissue_20221101115128_h0/simu_gwas_20221101115128',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.05/sim_tissue_20221101115128_h1/simu_gwas_20221101115128',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.05/sim_tissue_20221101115128_h2004/simu_gwas_20221101115128',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.05/sim_tissue_20221101115128_h20407/simu_gwas_20221101115128',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.05/sim_tissue_20221101115128_h20709/simu_gwas_20221101115128',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.2_hsq0.1/output_20000_20000_20221101115152/generated_20221101115152.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.1/sim_tissue_20221101115152_h0/simu_gwas_20221101115152',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.1/sim_tissue_20221101115152_h1/simu_gwas_20221101115152',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.1/sim_tissue_20221101115152_h2004/simu_gwas_20221101115152',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.1/sim_tissue_20221101115152_h20407/simu_gwas_20221101115152',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.1/sim_tissue_20221101115152_h20709/simu_gwas_20221101115152',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.2_hsq0.2/output_20000_20000_20221101115209/generated_20221101115209.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.2/sim_tissue_20221101115209_h0/simu_gwas_20221101115209',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.2/sim_tissue_20221101115209_h1/simu_gwas_20221101115209',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.2/sim_tissue_20221101115209_h2004/simu_gwas_20221101115209',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.2/sim_tissue_20221101115209_h20407/simu_gwas_20221101115209',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.2_hsq0.2/sim_tissue_20221101115209_h20709/simu_gwas_20221101115209',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.3_hsq0.05/output_20000_20000_20221101115240/generated_20221101115240.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.05/sim_tissue_20221101115240_h0/simu_gwas_20221101115240',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.05/sim_tissue_20221101115240_h1/simu_gwas_20221101115240',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.05/sim_tissue_20221101115240_h2004/simu_gwas_20221101115240',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.05/sim_tissue_20221101115240_h20407/simu_gwas_20221101115240',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.05/sim_tissue_20221101115240_h20709/simu_gwas_20221101115240',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.3_hsq0.1/output_20000_20000_20221101115256/generated_20221101115256.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.1/sim_tissue_20221101115256_h0/simu_gwas_20221101115256',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.1/sim_tissue_20221101115256_h1/simu_gwas_20221101115256',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.1/sim_tissue_20221101115256_h2004/simu_gwas_20221101115256',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.1/sim_tissue_20221101115256_h20407/simu_gwas_20221101115256',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.1/sim_tissue_20221101115256_h20709/simu_gwas_20221101115256',
        },
        {
            'generated_file_path': '/Volumes/HD/biodata/colocalization-tools/simu/risk1.3_hsq0.2/output_20000_20000_20221101115310/generated_20221101115310.tsv',
            'h0': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.2/sim_tissue_20221101115310_h0/simu_gwas_20221101115310',
            'h1': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.2/sim_tissue_20221101115310_h1/simu_gwas_20221101115310',
            'h2004': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.2/sim_tissue_20221101115310_h2004/simu_gwas_20221101115310',
            'h20407': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.2/sim_tissue_20221101115310_h20407/simu_gwas_20221101115310',
            'h20709': '/Volumes/HD/biodata/colocalization-tools/processed/risk1.3_hsq0.2/sim_tissue_20221101115310_h20709/simu_gwas_20221101115310',
        }
    ]

    rpt_set = []
    all_tool_list = [tool for tool, _, _ in TOOL_SIG_COL_INFO]
    for sim_files in sim_file_set:
        _rpts = {}
        for _sim_model in SIM_MODEL_ORDER:
            for tool in os.listdir(sim_files[_sim_model]):
                if tool not in all_tool_list:
                    continue
                if _rpts.get(tool) is None:
                    _rpts[tool] = {}
                _rpts[tool][_sim_model] = None
                analyzed_dir = os.path.join(sim_files[_sim_model], tool, 'analyzed')
                if not os.path.exists(analyzed_dir):
                    continue
                tool_sim_model_rpts = [os.path.join(analyzed_dir, tool_sim_model_rpt) for
                                       tool_sim_model_rpt in
                                       os.listdir(analyzed_dir)]
                tool_sim_model_rpts.sort(reverse=True, key=os.path.getmtime)
                if len(tool_sim_model_rpts) == 0:
                    continue
                for tool_sim_model_rpt in tool_sim_model_rpts:
                    if _rpts[tool][_sim_model] is not None:
                        break
                    if not (tool_sim_model_rpt.endswith('.tsv') or tool_sim_model_rpt.endswith('.tsv.gz')):
                        continue
                    if tool == 'fastenloc' and (tool_sim_model_rpt.endswith('.sig.tsv') or tool_sim_model_rpt.endswith('.sig.tsv.gz')):
                        _rpts[tool][_sim_model] = tool_sim_model_rpt
                    elif tool == 'ecaviar':
                        _rpts[tool][_sim_model] = tool_sim_model_rpt
                    elif os.path.basename(tool_sim_model_rpt).startswith(tool) and tool != 'fastenloc':
                        _rpts[tool][_sim_model] = tool_sim_model_rpt
        rpt_set.append((sim_files['generated_file_path'], _rpts))
    # training and testing model
    input_df_list = []
    for _idx, _rpt in enumerate(rpt_set):
        _input_df = prepare_sim_train_data(generated_file_path=_rpt[0], rpts=_rpt[1])
        _input_df[GENE_ID_COL_NAME] = _input_df[GENE_ID_COL_NAME] + f'#{os.path.basename(_rpt[0])}'
        input_df_list.append(_input_df)
    if len(input_df_list) == 0:
        print('No input files, nothing to do')
        exit(1)
    train_input_df = pd.concat(input_df_list)
    train_input_df.to_csv('training_input.tsv', sep='\t', header=True, index=False)
    train_model_xgb(train_input_df, 'test_set_result.tsv')
    train_model_xgb_rf(train_input_df, 'xgb_rf_test_set_result.tsv')
    train_svm_model(train_input_df, 'svm_test_set_result.tsv')
    train_rf_model(train_input_df, 'rf_test_set_result.tsv')
    train_lr_model(train_input_df, 'lr_test_set_result.tsv')
    train_pu_model(train_input_df, 'pu_test_set_result.tsv')
