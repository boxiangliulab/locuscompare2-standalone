import concurrent
import logging
import os
import sys
import traceback
import uuid
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

import pandas as pd

import common.config
import common.constants
from analyze.api import run_tools_api
from common import coloc_utils as utils, global_data_process as gdp
from figures import report_data_processor as redp
from ranking import scoring as sc


# def __before_run_jlim_tools_check(global_config):
#     logging.info(f'start check jlim')
#     utils.check_file_or_path_exist(global_config['input']['reference_genotype_panel_ld_ref_file'])


def __before_run_fastenloc_tools_check(global_config):
    logging.info(f'start check fastenloc')
    utils.check_file_or_path_exist(global_config['input']['ld_block_loci_file'])
    utils.check_file_or_path_exist(global_config['input']['eqtl_finemapping_file'])


def __before_run_coloc_tools_check(global_config):
    logging.info(f'start check coloc')


def __before_run_predixcan_tools_check(global_config):
    logging.info(f'start check predixcan')
    _model_db_path, _prediction_snp_covariance_path = utils.get_predixcan_ref_files(global_config)
    logging.info(f'Predixcan model db path: {_model_db_path}')
    logging.info(f'Predixcan SNP covariance path: {_prediction_snp_covariance_path}')
    utils.check_file_or_path_exist(_model_db_path)
    utils.check_file_or_path_exist(_prediction_snp_covariance_path)


def __before_run_smr_tools_check(global_config):
    logging.info(f'start check smr')
    utils.check_file_or_path_exist(global_config['input']['genecode'])
    __check_vcf(global_config)


def __before_run_ecaviar_tools_check(global_config):
    # no ref config to check
    logging.info(f'start check eCAVIAR')


def __before_run_twas_tools_check(global_config):
    logging.info(f'start check FUSION')
    twas_pos_path = utils.get_twas_ref_files(global_config)
    logging.info(f'TWAS pos file path: {twas_pos_path}')
    utils.check_file_or_path_exist(twas_pos_path)
    __check_vcf(global_config)


def __check_vcf(global_config):
    population = global_config.get('population', 'EUR').upper()
    ref_vcf_dir = global_config['input']['vcf']
    for chromosome in range(1, 23):
        input_vcf = os.path.join(ref_vcf_dir, population, f'chr{chromosome}.vcf.gz')
        utils.check_file_or_path_exist(input_vcf, False)


tools_func_map = {
    # 'jlim': {
    #     'check_fun': __before_run_jlim_tools_check,
    #     'run_fun': run_tools_api.__preprocess_and_run_jlim,
    # },
    'fastenloc': {
        'check_fun': __before_run_fastenloc_tools_check,
        'run_fun': run_tools_api.__preprocess_and_run_fastenloc,
    },
    'coloc': {
        'check_fun': __before_run_coloc_tools_check,
        'run_fun': run_tools_api.__preprocess_and_run_coloc,
    },
    'predixcan': {
        'check_fun': __before_run_predixcan_tools_check,
        'run_fun': run_tools_api.__preprocess_and_run_predixcan,
    },
    'smr': {
        'check_fun': __before_run_smr_tools_check,
        'run_fun': run_tools_api.__preprocess_and_run_smr,
    },
    'ecaviar': {
        'check_fun': __before_run_ecaviar_tools_check,
        'run_fun': run_tools_api.__preprocess_and_run_ecaviar,
    },
    'fusion': {
        'check_fun': __before_run_twas_tools_check,
        'run_fun': run_tools_api.__preprocess_and_run_twas,
    }
}


# def run(config_file=None, tools_list=None, log_file=None, parallel=False, tools_config=None, no_report=False):
def run(config_file=None, log_file=None, parallel=False, tools_config=None, no_report=False):

    # if tools_list is None:
    #     tools_list = ['all']
    cfg_list = []
    # retrieve config and study info
    if config_file is None:
        # default parameters
        study = common.constants.default_study
        cfg_list.append(common.constants.default_config)
    elif os.path.isdir(config_file):
        # specify study
        if config_file.endswith(os.sep):
            config_file = config_file.rstrip(os.sep)
        study = os.path.basename(config_file)
        for dir_path, _, file_names in os.walk(config_file):
            for cfg in file_names:
                if cfg.startswith('.'):
                    continue
                if not cfg.endswith('.yml') and not cfg.endswith('.yaml'):
                    continue
                cfg_list.append(os.path.join(dir_path, cfg))
    else:
        study = common.constants.default_study
        cfg_list.append(config_file)
    if len(cfg_list) == 0:
        logging.error(f'No config files to run.')
        return
        # loop to run all configs
    if log_file is None:
        log_file = f'{uuid.uuid4().hex}.log'

    report_list = []
    single_cfg_ensemble_result_ls = []
    numoftissues = len(cfg_list)
    currenttissuenum = 0
    for cfg in cfg_list:
        currenttissuenum = currenttissuenum + 1
        config_holder = common.config.ConfigHolder(single_config_file=cfg, 
                                                   study=study, parallel=parallel,
                                                   tools_config_file=tools_config)
        __init_logger(os.path.join(config_holder.study_dir, f'{log_file}'))
        # __run_single_cfg(tools_list, config_holder, report_list, parallel, study)
        single_cfg_ensemble_result = __run_single_cfg(config_holder, report_list, 
                                                      parallel, study, currenttissuenum, numoftissues)
        calculated_schedule = int(80/numoftissues * (currenttissuenum - 1))
        if os.path.exists('/process/'):
            with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
                schedule.write(str(calculated_schedule))
        else:
            with open(f"{os.path.join(config_holder.rank_dir, 'process_schedule.log')}", 'w') as schedule:
                schedule.write(str(calculated_schedule))
        schedule.close()
        if os.path.exists(single_cfg_ensemble_result): 
            if os.path.getsize(single_cfg_ensemble_result):
                single_cfg_ensemble_result_ls.append(single_cfg_ensemble_result)
        if len(single_cfg_ensemble_result_ls) > 1:
            union_result_gene(single_cfg_ensemble_result_ls)
        try:
            utils.cleanup_output(config_holder.tool_parent_dir)
        except:
            logging.warning(f'failed to clean {config_holder.tool_parent_dir}')
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(95))
    else:
        with open(f"{os.path.join(config_holder.rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(95))
    schedule.close()

    if len(report_list) == 0:
        logging.warning(f'No results of specified tools found')
        return
    if no_report is False:
        redp.report_data_process(report_list)


def __init_logger(logfile):
    log_format = '%(asctime)s::%(levelname)s::%(name)s::%(filename)s::%(lineno)d::%(message)s'
    stdout_hd = logging.StreamHandler(sys.stdout)
    stdout_hd.setLevel(logging.INFO)
    logging.root.handlers = []
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.FileHandler(logfile),
            stdout_hd,
        ]
    )


def get_tools_path(dir_path, tool_name):
    for path_name in os.listdir(os.path.join(dir_path, tool_name, 'analyzed')):
        if os.path.isfile(os.path.join(dir_path, tool_name, 'analyzed', path_name)) and (
                path_name.startswith(tool_name) or path_name.startswith('report')):
            return os.path.join(dir_path, tool_name, 'analyzed', path_name)

def union_result_gene(tissues_ensemble_result_ls):
    # 把每个tissue的基因统一
    all_tissue_df = pd.DataFrame()
    for single_cfg_ensemble_result in tissues_ensemble_result_ls:
        tmp = pd.read_csv(single_cfg_ensemble_result,sep='\t',usecols=['gene_id'])
        all_tissue_df = pd.concat([all_tissue_df,tmp],axis=0)
    all_tissue_df = all_tissue_df.dropna()
    all_gene = pd.DataFrame(list(set(all_tissue_df['gene_id'])))
    all_gene.columns = ['ensemble_gene_id']
    all_gene.index = all_gene['ensemble_gene_id']

    for single_cfg_ensemble_result in tissues_ensemble_result_ls:
        tmp = pd.read_csv(single_cfg_ensemble_result,sep='\t')
        tmp.index = tmp['gene_id']
        tmp = pd.concat([tmp,all_gene],axis=1)
        tmp['gene_id'] = tmp['ensemble_gene_id']
        tmp = tmp.iloc[:,:-1]
        tmp.to_csv(single_cfg_ensemble_result,sep='\t',index=False)
    pass


# def __run_single_cfg(tools_param_list, config_holder, report_list, parallel, study):
def __run_single_cfg(config_holder, report_list, parallel, study, currenttissuenum, numoftissues):
    start_time = datetime.now()
    tools_param_list = list(config_holder.global_config['tools'])
    logging.info(f'run tools_list: {tools_param_list}, start time: {start_time}')
    # check eqtl and gwas entry file exist
    utils.check_file_or_path_exist(config_holder.global_config['working_dir'])
    utils.check_file_or_path_exist(config_holder.global_config['input']['gwas']['file'])
    utils.check_file_or_path_exist(config_holder.global_config['input']['eqtl']['file'])
    # utils.check_file_or_path_exist(global_config['input']['vcf'])

    # check tool require file exist
    print(type(tools_param_list))
    actually_tools_list = []
    if 'all' in tools_param_list:
        actually_tools_list = tools_func_map.keys()
    else:
        for tool in tools_param_list:
            print(tool)
            if tools_func_map.get(tool):
                actually_tools_list.append(tool)
            else:
                logging.error(f'The {tool} tool is not recognized')
    print("*****actually_tools_list*****")
    print(actually_tools_list)

    smr_schedual = False
    ecaviar_schedual = False
    coloc_schedual = False
    fastenloc_schedual = False
    fusion_schedual = False
    predixcan_schedual = False

    for tool in actually_tools_list:
        if tool == 'smr':
            smr_schedual = True
            break
        if tool == 'ecaviar':
            ecaviar_schedual = True
            break
        if tool == 'coloc':
            coloc_schedual = True
            break
        if tool == 'fastenloc':
            fastenloc_schedual = True
            break
        if tool == 'predixcan':
            predixcan_schedual = True
            break
        if tool == 'fusion':
            fusion_schedual = True
            break

        
    config_schedual = {'smr': smr_schedual,
                        'ecaviar': ecaviar_schedual,
                        'coloc': coloc_schedual,
                        'fastenloc': fastenloc_schedual,
                        'fusion': fusion_schedual,
                        'predixcan': predixcan_schedual}

    # check tool require file exist
    for tool in actually_tools_list:
        tools_func_map[tool]['check_fun'](config_holder.global_config)

    # preprocess gwas,eqtl,vcf
    processor = gdp.Processor(config_holder)
    if not utils.file_exists(processor.gwas_preprocessed_file):
        processor.preprocess_gwas()
    gwas_sig_df = pd.read_table(config_holder.gwas_filter_file, nrows=2)
    if gwas_sig_df.shape[0] == 0:
        gwas_file = config_holder.global_config['input']['gwas']['file']
        logging.warning(f'No significant records found in GWAS file {gwas_file} '
                        f'by threshold {config_holder.gwas_p_threshold}')
        return
    if not utils.file_exists(processor.eqtl_output_report):
        processor.preprocess_eqtl()
    eqtl_sig_df = pd.read_table(config_holder.eqtl_output_report, nrows=2)
    if eqtl_sig_df.shape[0] == 0:
        eql_file = config_holder.global_config['input']['eqtl']['file']
        logging.warning(f'No significant records found in eQTL file {eql_file} '
                        f'by threshold {config_holder.eqtl_p_threshold}')
        return
    # check preprocess file is exist
    utils.check_path_exist_and_has_size(processor.eqtl_output_report)
    utils.check_path_exist_and_has_size(processor.gwas_preprocessed_file)
    utils.check_file_or_path_exist(processor.gwas_cluster_output_dir)
    results = dict()
    rank_output_file = os.path.join(processor.rank_dir,
                                    f'ensemble_ranking_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv')
    if parallel:
        logging.info("Parallel running {}".format(actually_tools_list))
        with ProcessPoolExecutor(max_workers=len(actually_tools_list)) as executor:
            tool_futures = {}
            for tool in actually_tools_list:
                tool_futures[executor.submit(tools_func_map[tool]['run_fun'], processor, currenttissuenum, numoftissues, whether_schedual = config_schedual[tool])] = tool
            exceptions = []
            for future in concurrent.futures.as_completed(tool_futures.keys()):
                current_tool = tool_futures[future]
                try:
                    results[current_tool] = future.result()
                    logging.info(f'Tool {current_tool} completed successfully!')
                except Exception as exc:
                    exceptions.append(exc)
                    logging.info(f'Tool {current_tool} failed with error: {exc}!')
            if exceptions:
                for ex in exceptions:
                    logging.error("".join(traceback.TracebackException.from_exception(ex).format()))
                raise Exception(exceptions)
        for tool in results.keys():
            report_list.append(
                {'trait': config_holder.gwas_trait, 'tool_name': tool, 'study1': study,
                 'tissue': config_holder.eqtl_tissue, 'report_path': results[tool],
                 'cfg_pro': processor, 'rank_output_file': rank_output_file})
    else:
        for tool in actually_tools_list:
            logging.info(f'Tool {tool} start running')
            try:
                report_file_path = tools_func_map[tool]['run_fun'](processor, currenttissuenum, numoftissues, whether_schedual = config_schedual[tool])
                logging.info(f'Tool {tool} completed successfully!')
            except Exception as error:
                logging.info(f'Tool {tool} failed with error: {error}!')
                raise error
            results[tool] = report_file_path
            report_list.append({'trait': config_holder.gwas_trait, 'tool_name': tool, 'study1': study,
                                'tissue': config_holder.eqtl_tissue, 'report_path': report_file_path,
                                'cfg_pro': processor, 'rank_output_file': rank_output_file})
    sc.run_ranking(rpt_obj=results, output_file_path=rank_output_file,
                   sample_size=processor.global_config['input']['gwas']['sample_size'])
    logging.info(f'coloctools complete at: {datetime.now()},duration: {datetime.now() - start_time}')
    logging.info(f'tissue: {config_holder.eqtl_tissue},trait: {config_holder.gwas_trait} coloctools complete')
    return rank_output_file


if __name__ == '__main__':
    start_time = datetime.now()
    logging.info(f'start run all coloctools, start time: {start_time}')
    parse_args = utils.parse_parameters()
    # run(parse_args.config_file, parse_args.tools_list, parse_args.log_file, parse_args.parallel,
    #     parse_args.tools_config, parse_args.no_report)
    # run(parse_args.config_file, parse_args.log_file, parse_args.parallel,
    run(parse_args.config_file, parse_args.log_file, True,
        parse_args.tools_config, parse_args.no_report)
    logging.info(f'all coloctools complete at: {datetime.now()},duration: {datetime.now() - start_time}')
