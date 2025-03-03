import asyncio
import os
from pathlib import Path
import sys

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
print(sys.path)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print(sys.path)

import eQTL.coloc.run_coloc as rc
import eQTL.fastenloc.run_fastenloc as rf
import eQTL.fastenloc.gwas_data_processor as fgdp
from eQTL.predixcan import run_predixcan as rp, gwas_data_processor as pgdp
from eQTL.smr import run_smr as rs, qtl_data_processor as sedp
from eQTL.ecaviar import run_ecaviar as run_e, data_processor as edp
from eQTL.twas import run_twas as rt
from common import global_data_process as gdp
from common import utils
import logging



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
    gwas_threshold = glob_processor.config_holder.gwas_p_threshold
    qtl_threshold = glob_processor.config_holder.qtl_p_threshold
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
                     min_matching_number = glob_processor.min_matching_number,
                     gwas_threshold = gwas_threshold,
                     qtl_threshold = qtl_threshold)


######################
###  2. fastenloc  ###
######################

def __preprocess_and_run_fastenloc(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_fastenloc: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    _working_dir = os.path.join(glob_processor.tool_parent_dir, rf.Fastenloc.COLOC_TOOL_NAME)
    fgdp_obj = fgdp.FastenlocGwasProcessor()
    # fastenloc_gwas_result, gwas_snp_count = fgdp_obj.prepare_gwas_data(\
    #                 working_dir=_working_dir,
    #                 gwas_preprocessed_file=glob_processor.gwas_preprocessed_file,
    #                 gwas_col_dict=glob_processor.gwas_col_dict,
    #                 ld_block_loci_file=glob_processor.global_config['input']['ld_block_loci_file'], 
    #                 rank_dir = glob_processor.rank_dir,
    #                 current_analysis_order=current_analysis_order, total_numof_analyses=total_numof_analyses, 
    #                 whether_schedual=whether_schedual)

    # rf_obj = rf.Fastenloc()
    # return rf_obj.run(qtl_tissue=glob_processor.qtl_tissue,
    #                   working_dir=_working_dir,
    #                   qtl_finemapping_file=glob_processor.global_config['input']['qtl_finemapping_file'],
    #                   qtl_output_report=glob_processor.qtl_output_report,
    #                   output_torus_output_file=fastenloc_gwas_result,
    #                   gwas_snp_count=gwas_snp_count,
    #                   tools_config_file=glob_processor.tools_config_file)

    gwas_threshold = glob_processor.config_holder.gwas_p_threshold
    qtl_threshold = glob_processor.config_holder.qtl_p_threshold

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
                     min_matching_number = glob_processor.min_matching_number,
                     gwas_threshold = gwas_threshold,
                     qtl_threshold = qtl_threshold)

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

    gwas_threshold = glob_processor.config_holder.gwas_p_threshold
    qtl_threshold = glob_processor.config_holder.qtl_p_threshold

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
                            min_matching_number = glob_processor.min_matching_number,
                            gwas_threshold = gwas_threshold,
                            qtl_threshold = qtl_threshold)

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
    print(f"glob_processor.ref_vcf_dir: {glob_processor.ref_vcf_dir}")
    # prepare gwas file finished
    pop = glob_processor.global_config.get('population', 'EUR').upper()

    gwas_threshold = glob_processor.config_holder.gwas_p_threshold
    qtl_threshold = glob_processor.config_holder.qtl_p_threshold

    smr_qtl_processor = sedp.SmrQTLProcessor()
    smr_qtl_result = smr_qtl_processor.prepare_ld_ref(working_dir = _working_dir,
                                                        var_id_col_name = gdp.Processor.VAR_ID_COL_NAME,
                                                        gwas_filter_file = glob_processor.gwas_filter_file,
                                                        gwas_col_dict = glob_processor.gwas_col_dict,
                                                        qtl_output_report = glob_processor.qtl_output_report,
                                                        qtl_grouped_dir = glob_processor.qtl_grouped_dir,
                                                        qtl_col_dict = glob_processor.qtl_col_dict,
                                                        ref_vcf_dir = glob_processor.ref_vcf_dir,
                                                        population = pop,
                                                        rank_dir = glob_processor.rank_dir,
                                                        current_analysis_order = current_analysis_order,
                                                        total_numof_analyses = total_numof_analyses,
                                                        whether_schedual = whether_schedual,
                                                        min_matching_number = glob_processor.min_matching_number,
                                                        gwas_threshold = gwas_threshold,
                                                        qtl_threshold = qtl_threshold)
    
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

    return smr.run(working_dir = _working_dir,
                   var_id_col_name = gdp.Processor.VAR_ID_COL_NAME,
                   genecode_file = _genecode_file,
                   gwas_filter_file = glob_processor.gwas_filter_file,
                   gwas_chrom_group_dir = glob_processor.gwas_output_dir,
                   gwas_col_dict = glob_processor.gwas_col_dict,
                   gwas_sample_size = _gwas_sample_size,
                   qtl_output_report = glob_processor.qtl_output_report,
                   qtl_grouped_dir = glob_processor.qtl_grouped_dir,
                   qtl_col_dict = glob_processor.qtl_col_dict,
                   subset_vcf_dir = _subset_vcf_dir,
                   subset_vcf_file_pattern = smr_qtl_result[1],
                   ld_ref_dir = _ld_ref_dir,
                   ld_ref_file_pattern = smr_qtl_result[3],
                   tools_config_file = glob_processor.tools_config_file,
                   min_matching_number = glob_processor.min_matching_number)

######################
###  5. predixcan  ###
######################

def __preprocess_and_run_predixcan(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_predixcan: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    _working_dir = os.path.join(glob_processor.tool_parent_dir, pgdp.PredixcanGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    gwas_data_processor = pgdp.PredixcanGwasProcessor()
    ## whole gwas
    _gwas_processed_file = gwas_data_processor.prepare(_working_dir,
                                                       glob_processor.gwas_preprocessed_file,
                                                       gdp.Processor.VAR_ID_COL_NAME,
                                                       glob_processor.gwas_col_dict,
                                                       rank_dir = glob_processor.rank_dir,
                                                       current_analysis_order = current_analysis_order, 
                                                       total_numof_analyses = total_numof_analyses, 
                                                       whether_schedual = whether_schedual)
    
    ## gwas loci
    # _gwas_processed_file = gwas_data_processor.prepare(_working_dir,
    #                                                    glob_processor.gwas_cluster_output_dir,
    #                                                    gdp.Processor.VAR_ID_COL_NAME,
    #                                                    glob_processor.gwas_col_dict,
    #                                                    rank_dir = glob_processor.rank_dir,
    #                                                    current_analysis_order = current_analysis_order, 
    #                                                    total_numof_analyses = total_numof_analyses, 
    #                                                    whether_schedual = whether_schedual)
    
    
    # prepare gwas file finished
    # if not os.path.exists(_gwas_processed_file) or os.path.getsize(_gwas_processed_file) <= 0:
    #     logging.warning(f'Dependent files not found, did you run gwas_data_processor?')
    #     return None
    if len(_gwas_processed_file) == 0:
        return None
    _model_db_path, _prediction_snp_covariance_path = utils.get_predixcan_ref_files(
        glob_processor.config_holder.global_config)
    predixcan = rp.Predixcan()
    gwas_sample_size = int(glob_processor.global_config['input']['gwas']['sample_size'])
    return predixcan.run(_working_dir,
                         _model_db_path,
                         _prediction_snp_covariance_path,
                         _gwas_processed_file,
                         glob_processor.gwas_col_dict,
                         glob_processor.qtl_output_report, gwas_sample_size=gwas_sample_size)

######################
###    6. fusion   ###
######################

def __preprocess_and_run_twas(glob_processor, current_analysis_order, total_numof_analyses, whether_schedual):
    print(f"__preprocess_and_run_twas: {glob_processor.global_config['input']['qtl']['qtl_type']}")

    ## chrom
    _working_dir = os.path.join(glob_processor.tool_parent_dir, rt.TWAS.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    if len(os.listdir(glob_processor.gwas_output_dir)) == 0:
        raise ValueError(f'Dependant GWAS files not found, did you run global_data_process?')
    _weight_pos_file = utils.get_twas_ref_files(glob_processor.global_config)
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    twas = rt.TWAS()

    return twas.run(_working_dir,
                    _weight_pos_file,
                    glob_processor.gwas_output_dir,
                    glob_processor.gwas_col_dict,
                    glob_processor.ref_vcf_dir,
                    pop,
                    parallel=False,
                    tools_config_file=glob_processor.tools_config_file)


    ## loci 
    # _working_dir = os.path.join(glob_processor.tool_parent_dir, rt.TWAS.COLOC_TOOL_NAME)
    # Path(_working_dir).mkdir(exist_ok=True, parents=True)
    # if len(os.listdir(glob_processor.gwas_cluster_output_dir)) == 0:
    #     raise ValueError(f'Dependant GWAS files not found, did you run global_data_process?')
    # _weight_pos_file = utils.get_twas_ref_files(glob_processor.global_config)
    # pop = glob_processor.global_config.get('population', 'EUR').upper()
    # twas = rt.TWAS()

    # return twas.run(_working_dir,
    #                 _weight_pos_file,
    #                 glob_processor.gwas_cluster_output_dir,
    #                 glob_processor.gwas_col_dict,
    #                 glob_processor.ref_vcf_dir,
    #                 pop,
    #                 parallel=False,
    #                 tools_config_file=glob_processor.tools_config_file)