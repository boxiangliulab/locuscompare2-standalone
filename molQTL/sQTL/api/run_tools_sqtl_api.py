import asyncio
import os
from pathlib import Path
import sQTL.coloc.run_coloc as rc
import sQTL.fastenloc.run_fastenloc as rf
import sQTL.fastenloc.gwas_data_processor as fgdp
from sQTL.predixcan import run_predixcan as rp, gwas_data_processor as pgdp
from sQTL.smr import run_smr as rs, qtl_data_processor as sedp
from sQTL.ecaviar import run_ecaviar as run_e, data_processor as edp
from sQTL.twas import run_twas as rt
from common import global_data_process as gdp
from common import utils
import logging
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print(f"sqtl sys path {sys.path}")

######################
###    1. coloc    ###
######################

def __preprocess_and_run_coloc(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_coloc: {glob_processor.global_config['input']['qtl']['qtl_type']}")
    _working_dir = os.path.join(glob_processor.tool_parent_dir, rc.Coloc.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _qtl_sample_size = glob_processor.global_config['input']['qtl'].get('sample_size', 948)
    _gwas_type = glob_processor.global_config['input']['gwas'].get('type', 'cc')
    _qtl_type = glob_processor.global_config['input']['qtl'].get('type', 'quant')
    coloc = rc.Coloc()

    return coloc.run(_working_dir,
                     gdp.Processor.VAR_ID_COL_NAME,
                     glob_processor.gwas_cluster_output_dir,
                     glob_processor.gwas_cluster_summary,
                     glob_processor.gwas_col_dict,
                     _gwas_sample_size,
                     glob_processor.qtl_output_report,
                     glob_processor.qtl_grouped_dir,
                     glob_processor.qtl_col_dict,
                     _qtl_sample_size,
                     _gwas_type,
                     _qtl_type,
                     tools_config=glob_processor.tools_config_file,
                     parallel=glob_processor.config_holder.parallel,
                     rank_dir = glob_processor.rank_dir,
                     current_analysis_order = current_analysis_order, 
                     total_numof_analyses = total_numof_analyses, 
                     whether_schedual = whether_schedual,
                     min_matching_number = glob_processor.min_matching_number)


######################
###  2. fastenloc  ###
######################


def __preprocess_and_run_fastenloc(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_fastenloc: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    _working_dir = os.path.join(glob_processor.tool_parent_dir, rf.Fastenloc.COLOC_TOOL_NAME)
    fgdp_obj = fgdp.FastenlocGwasProcessor()

    output_summary_statistics_file = fgdp_obj.prepare_fastenlocinput_data(_working_dir,
                     gdp.Processor.VAR_ID_COL_NAME,
                     glob_processor.gwas_cluster_output_dir,
                     glob_processor.gwas_cluster_summary,
                     glob_processor.gwas_col_dict,
                     glob_processor.qtl_output_report,
                     glob_processor.qtl_grouped_dir,
                     glob_processor.qtl_col_dict,
                     parallel=glob_processor.config_holder.parallel,
                     rank_dir = glob_processor.rank_dir,
                     current_analysis_order = current_analysis_order, 
                     total_numof_analyses = total_numof_analyses, 
                     whether_schedual = whether_schedual,
                     min_matching_number = glob_processor.min_matching_number)

    rf_obj = rf.Fastenloc()
    return rf_obj.run(biological_context=glob_processor.biological_context,
                      working_dir=_working_dir,
                      qtl_output_report = glob_processor.qtl_output_report,
                      summary_statistics_file=output_summary_statistics_file,
                      tools_config_file=glob_processor.tools_config_file)


######################
###   3. ecaviar   ###
######################

def __preprocess_and_run_ecaviar(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_ecaviar: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    _working_dir = os.path.join(glob_processor.tool_parent_dir, run_e.ECaviar.ECAVIAR_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _qtl_sample_size = glob_processor.global_config['input']['qtl'].get('sample_size', 948)

    ecaviar_data_processor = edp.ECaviarDataProcessor()
    preproc_rst_dir = ecaviar_data_processor.prepare(\
                            working_dir=_working_dir,
                            gwas_cluster_dir=glob_processor.gwas_cluster_output_dir,
                            gwas_cluster_summary=glob_processor.gwas_cluster_summary,
                            qtl_group_dir=glob_processor.qtl_grouped_dir,
                            qtl_report=glob_processor.qtl_output_report,
                            ref_vcf_dir=glob_processor.ref_vcf_dir,
                            gwas_col_dict=glob_processor.gwas_col_dict,
                            qtl_col_dict=glob_processor.qtl_col_dict,
                            population=glob_processor.global_config.get('population',
                                                                        'EUR').upper(),
                            gwas_sample_size=_gwas_sample_size,
                            qtl_sample_size=_qtl_sample_size,
                            var_id_col_name=gdp.Processor.VAR_ID_COL_NAME,
                            rank_dir = glob_processor.rank_dir,
                            current_analysis_order = current_analysis_order, 
                            total_numof_analyses = total_numof_analyses, 
                            whether_schedual = whether_schedual,
                            min_matching_number = glob_processor.min_matching_number)

    ecaviar = run_e.ECaviar()
    return asyncio.run(ecaviar.run(working_dir=_working_dir,
                                   candidate_data_dir=preproc_rst_dir,
                                   parallel=glob_processor.config_holder.parallel,
                                   tools_config=glob_processor.tools_config_file,
                                   rank_dir = glob_processor.rank_dir,
                                   current_analysis_order = current_analysis_order, 
                                   total_numof_analyses = total_numof_analyses, 
                                   whether_schedual = whether_schedual))


######################
###     4. smr     ###
######################

def __preprocess_and_run_smr(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_smr: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    _working_dir = os.path.join(glob_processor.tool_parent_dir, rs.Smr.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    # prepare gwas file finished
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    smr_qtl_processor = sedp.SmrQTLProcessor()
    smr_qtl_result = smr_qtl_processor.prepare_ld_ref(_working_dir,
                                                        gdp.Processor.VAR_ID_COL_NAME,
                                                        glob_processor.gwas_filter_file,
                                                        glob_processor.gwas_col_dict,
                                                        glob_processor.qtl_output_report,
                                                        glob_processor.qtl_grouped_dir,
                                                        glob_processor.qtl_col_dict,
                                                        glob_processor.config_holder.qtl_p_threshold,
                                                        glob_processor.ref_vcf_dir,
                                                        pop,
                                                        glob_processor.rank_dir,
                                                        current_analysis_order,
                                                        total_numof_analyses,
                                                        whether_schedual,
                                                        min_matching_number = glob_processor.min_matching_number)
    # prepare ldref file finished
    if len(os.listdir(glob_processor.gwas_output_dir)) == 0:
        logging.warning('Dependent files not found, did you run gwas_data_processor?')
        return None
    _subset_vcf_dir = smr_qtl_result[0]
    _ld_ref_dir = smr_qtl_result[2]
    if len(os.listdir(_ld_ref_dir)) == 0 or len(os.listdir(_subset_vcf_dir)) == 0:
        logging.warning('Dependant files not found, did you run qtl_data_processor?')
        return None

    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _qtl_sample_size = glob_processor.global_config['input']['qtl'].get('sample_size', 948)
    _genecode_file = glob_processor.global_config['input']['genecode']
    smr = rs.Smr()

    return smr.run(_working_dir,
                   gdp.Processor.VAR_ID_COL_NAME,
                   _genecode_file,
                   glob_processor.gwas_filter_file,
                   glob_processor.gwas_output_dir,
                   glob_processor.gwas_col_dict,
                   _gwas_sample_size,
                   glob_processor.qtl_output_report,
                   glob_processor.qtl_grouped_dir,
                   glob_processor.qtl_col_dict,
                   glob_processor.config_holder.qtl_p_threshold,
                   _subset_vcf_dir,
                   smr_qtl_result[1],
                   _ld_ref_dir,
                   smr_qtl_result[3],
                   tools_config_file=glob_processor.tools_config_file,
                   min_matching_number = glob_processor.min_matching_number)


######################
###  5. predixcan  ###
######################

def __preprocess_and_run_predixcan(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_predixcan: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    _working_dir = os.path.join(glob_processor.tool_parent_dir, pgdp.PredixcanGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    gwas_data_processor = pgdp.PredixcanGwasProcessor()
    _gwas_processed_file = gwas_data_processor.prepare(_working_dir,
                                                       glob_processor.gwas_preprocessed_file,
                                                       gdp.Processor.VAR_ID_COL_NAME,
                                                       glob_processor.gwas_col_dict,
                                                       rank_dir = glob_processor.rank_dir,
                                                       current_analysis_order = current_analysis_order, 
                                                       total_numof_analyses = total_numof_analyses, 
                                                       whether_schedual = whether_schedual)
    
    # prepare gwas file finished
    if not os.path.exists(_gwas_processed_file) or os.path.getsize(_gwas_processed_file) <= 0:
        logging.warning(f'Dependent files not found, did you run gwas_data_processor?')
        return None
    _model_db_path, _prediction_snp_covariance_path = utils.get_predixcan_ref_files(
        glob_processor.config_holder.global_config)
    predixcan = rp.Predixcan()
    return predixcan.run(_working_dir,
                         _model_db_path,
                         _prediction_snp_covariance_path,
                         _gwas_processed_file,
                         glob_processor.gwas_col_dict,
                         glob_processor.qtl_output_report)


######################
###    6. fusion   ###
######################

def __preprocess_and_run_twas(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_twas: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    _working_dir = os.path.join(glob_processor.tool_parent_dir, rt.TWAS.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    if len(os.listdir(glob_processor.gwas_output_dir)) == 0:
        raise ValueError(f'Dependant GWAS files not found, did you run global_data_process?')
    _weight_pos_file = utils.get_twas_ref_files(glob_processor.global_config)
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    twas = rt.TWAS()
    # parallel = glob_processor.config_holder.parallel
    # if parallel:
    #     memory_size = psutil.virtual_memory().total / 1024 / 1024 / 1024
    #     # baseline is 16GB, every TWAS worker use up to 5GB memory
    #     worker_num = int((memory_size - 16) // 5)
    #     if worker_num < 2:
    #         worker_num = 1
    #         parallel = False
    #     elif worker_num > 22:
    #         worker_num = 22
    # else:
    #     worker_num = 1
    return twas.run(_working_dir,
                    _weight_pos_file,
                    glob_processor.gwas_output_dir,
                    glob_processor.gwas_col_dict,
                    glob_processor.ref_vcf_dir,
                    pop,
                    parallel=False,
                    tools_config_file=glob_processor.tools_config_file)
