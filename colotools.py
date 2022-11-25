import os
from multiprocessing import Pool

from common import coloc_utils as utils, global_data_process as gdp
from analyze.api import run_tools_api
from figures import report_data_processor as redp
from datetime import datetime
import common.config
import common.constants
import uuid
import logging
import sys


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
    utils.check_file_or_path_exist(_model_db_path)
    utils.check_file_or_path_exist(_prediction_snp_covariance_path)


def __before_run_smr_tools_check(global_config):
    logging.info(f'start check smr')
    utils.check_file_or_path_exist(global_config['input']['genecode'])
    __check_vcf(global_config)


def __before_run_ecaviar_tools_check(global_config):
    # no ref config to check
    logging.info(f'start check eCAVIAR')


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
    }
}


def run(config_file=None, tools_list=None, log_file=None, parallel=False):
    if tools_list is None:
        tools_list = ['all']
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
    for cfg in cfg_list:
        config_holder = common.config.ConfigHolder(single_config_file=cfg, study=study, parallel=parallel)
        __init_logger(os.path.join(config_holder.study_dir, f'{log_file}'))
        __run_single_cfg(tools_list, config_holder, report_list, parallel, study)
        # todo create a process for the clean up job
        # try:
        #     utils.cleanup_output(config_holder.tool_parent_dir)
        # except:
        #     logging.warning(f'failed to clean {config_holder.tool_parent_dir}')
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


def __run_single_cfg(tools_param_list, config_holder, report_list, parallel, study):
    start_time = datetime.now()
    logging.info(f'run tools_list: {tools_param_list}, start time: {start_time}')
    # check eqtl and gwas entry file exist
    utils.check_file_or_path_exist(config_holder.global_config['working_dir'])
    utils.check_file_or_path_exist(config_holder.global_config['input']['gwas']['file'])
    utils.check_file_or_path_exist(config_holder.global_config['input']['eqtl']['file'])
    # utils.check_file_or_path_exist(global_config['input']['vcf'])

    # check tool require file exist
    actually_tools_list = []
    if 'all' in tools_param_list:
        actually_tools_list = tools_func_map.keys()
    else:
        for tool in tools_param_list:
            if tools_func_map.get(tool):
                actually_tools_list.append(tool)
            else:
                logging.error(f'The {tool} tool is not recognized')

    # check tool require file exist
    for tool in actually_tools_list:
        tools_func_map[tool]['check_fun'](config_holder.global_config)

    # preprocess gwas,eqtl,vcf
    processor = gdp.Processor(config_holder)
    if not utils.file_exists(processor.gwas_preprocessed_file):
        processor.preprocess_gwas()
    if not utils.file_exists(processor.eqtl_output_report):
        processor.preprocess_eqtl()

    # check preprocess file is exist
    utils.check_path_exist_and_has_size(processor.eqtl_output_report)
    utils.check_path_exist_and_has_size(processor.gwas_preprocessed_file)
    utils.check_file_or_path_exist(processor.gwas_cluster_output_dir)

    if parallel:
        logging.info("Parallel run {}".format(actually_tools_list))
        p = Pool(len(actually_tools_list))

        results = dict()
        for tool in actually_tools_list:
            logging.info("Run tool {}".format(tool))
            results[tool] = p.apply_async(tools_func_map[tool]['run_fun'], args=(processor,))

        p.close()
        p.join()
        for tool in results.keys():
            report_list.append(
                {'trait': config_holder.gwas_trait, 'tool_name': tool, 'study1': study,
                 'tissue': config_holder.eqtl_tissue, 'report_path': results[tool].get(),
                 'cfg_pro': processor})
    else:
        for tool in actually_tools_list:
            logging.info("Run tool {}".format(tool))
            # report_file_path = get_tools_path(processor.tool_parent_dir, tool)
            report_file_path = tools_func_map[tool]['run_fun'](processor)
            report_list.append({'trait': config_holder.gwas_trait, 'tool_name': tool, 'study1': study,
                                'tissue': config_holder.eqtl_tissue, 'report_path': report_file_path,
                                'cfg_pro': processor})

    logging.info(f'coloctools complete at: {datetime.now()},duration: {datetime.now() - start_time}')


if __name__ == '__main__':
    parse_args = utils.parse_parameters()
    run(parse_args.config_file, parse_args.tools_list, parse_args.log_file, parse_args.parallel)
