import asyncio
import os
from pathlib import Path
import analyze.jlim.run_jlim as jt
import analyze.coloc.run_coloc as rc
import analyze.fastenloc.run_fastenloc as rf
import analyze.fastenloc.gwas_data_processor as fgdp
from analyze.predixcan import run_predixcan as rp, gwas_data_processor as pgdp
from analyze.smr import run_smr as rs, eqtl_data_processor as sedp
from analyze.ecaviar import run_ecaviar as run_e, data_processor as edp
from analyze.twas import run_twas as rt
from common import global_data_process as gdp
from common import coloc_utils as util
import logging


def __preprocess_and_run_jlim(global_processor):
    global_config = global_processor.global_config
    eqtl_file_path = f'{global_config["input"]["eqtl"]["file"]}'
    reference_genotype_panel_ld_ref_file = global_config['input']['reference_genotype_panel_ld_ref_file']
    eqtl_sample_size = global_config['input']['eqtl']['sample_size']
    _working_dir = os.path.join(global_processor.tool_parent_dir, jt.Jlim.COLOC_TOOL_NAME)
    jlim_obj = jt.Jlim()
    return jlim_obj.run(gwas_col_dict=global_processor.gwas_col_dict,
                        eqtl_col_dict=global_processor.eqtl_col_dict,
                        var_id_col_name=gdp.Processor.VAR_ID_COL_NAME,
                        eqtl_file_path=eqtl_file_path,
                        reference_genotype_panel_ld_ref_file=reference_genotype_panel_ld_ref_file,
                        eqtl_sample_size=eqtl_sample_size,
                        working_dir=_working_dir,
                        gwas_cluster_output_dir=global_processor.gwas_cluster_output_dir,
                        eqtl_output_report=global_processor.eqtl_output_report,
                        gwas_filter_file=global_processor.gwas_filter_file,
                        eqtl_output_dir=global_processor.eqtl_output_dir,
                        parallel=global_processor.config_holder.parallel)


def __preprocess_and_run_fastenloc(processor):
    _working_dir = os.path.join(processor.tool_parent_dir, rf.Fastenloc.COLOC_TOOL_NAME)
    fgdp_obj = fgdp.FastenlocGwasProcessor()
    fastenloc_gwas_result, gwas_snp_count = fgdp_obj.prepare_gwas_data(working_dir=_working_dir,
                                                       gwas_preprocessed_file=processor.gwas_preprocessed_file,
                                                       gwas_col_dict=processor.gwas_col_dict,
                                                       ld_block_loci_file=processor.global_config['input'][
                                                           'ld_block_loci_file'])

    rf_obj = rf.Fastenloc()
    return rf_obj.run(eqtl_tissue=processor.eqtl_tissue,
                      working_dir=_working_dir,
                      eqtl_finemapping_file=processor.global_config['input']['eqtl_finemapping_file'],
                      eqtl_output_report=processor.eqtl_output_report,
                      output_torus_output_file=fastenloc_gwas_result,
                      gwas_snp_count=gwas_snp_count,
                      tools_config_file=processor.tools_config_file)


def __preprocess_and_run_coloc(glob_processor):
    _working_dir = os.path.join(glob_processor.tool_parent_dir, rc.Coloc.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _eqtl_sample_size = glob_processor.global_config['input']['eqtl'].get('sample_size', 948)
    _gwas_type = glob_processor.global_config['input']['gwas'].get('type', 'cc')
    _eqtl_type = glob_processor.global_config['input']['eqtl'].get('type', 'quant')
    coloc = rc.Coloc()

    return coloc.run(_working_dir,
                     gdp.Processor.VAR_ID_COL_NAME,
                     glob_processor.gwas_cluster_output_dir,
                     glob_processor.gwas_cluster_summary,
                     glob_processor.gwas_col_dict,
                     _gwas_sample_size,
                     glob_processor.eqtl_output_report,
                     glob_processor.eqtl_output_dir,
                     glob_processor.eqtl_col_dict,
                     _eqtl_sample_size,
                     _gwas_type,
                     _eqtl_type,
                     tools_config=glob_processor.tools_config_file,
                     parallel=glob_processor.config_holder.parallel)


def __preprocess_and_run_predixcan(glob_processor):
    _working_dir = os.path.join(glob_processor.tool_parent_dir, pgdp.PredixcanGwasProcessor.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    gwas_data_processor = pgdp.PredixcanGwasProcessor()
    _gwas_processed_file = gwas_data_processor.prepare(_working_dir,
                                                       glob_processor.gwas_preprocessed_file,
                                                       gdp.Processor.VAR_ID_COL_NAME,
                                                       glob_processor.gwas_col_dict)
    # prepare gwas file finished
    if not os.path.exists(_gwas_processed_file) or os.path.getsize(_gwas_processed_file) <= 0:
        logging.warning(f'Dependent files not found, did you run gwas_data_processor?')
        return None
    _model_db_path, _prediction_snp_covariance_path = util.get_predixcan_ref_files(
        glob_processor.config_holder.global_config)
    predixcan = rp.Predixcan()
    return predixcan.run(_working_dir,
                         _model_db_path,
                         _prediction_snp_covariance_path,
                         _gwas_processed_file,
                         glob_processor.gwas_col_dict,
                         glob_processor.eqtl_output_report)


def __preprocess_and_run_smr(glob_processor):
    _working_dir = os.path.join(glob_processor.tool_parent_dir, rs.Smr.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    # prepare gwas file finished
    pop = glob_processor.global_config.get('population', 'EUR').upper()
    smr_eqtl_processor = sedp.SmrEqtlProcessor()
    smr_eqtl_result = smr_eqtl_processor.prepare_ld_ref(_working_dir,
                                                        gdp.Processor.VAR_ID_COL_NAME,
                                                        glob_processor.gwas_filter_file,
                                                        glob_processor.gwas_col_dict,
                                                        glob_processor.eqtl_output_report,
                                                        glob_processor.eqtl_output_dir,
                                                        glob_processor.eqtl_col_dict,
                                                        glob_processor.config_holder.eqtl_p_threshold,
                                                        glob_processor.ref_vcf_dir,
                                                        pop)
    # prepare ldref file finished
    if len(os.listdir(glob_processor.gwas_output_dir)) == 0:
        logging.warning('Dependent files not found, did you run gwas_data_processor?')
        return None
    _subset_vcf_dir = smr_eqtl_result[0]
    _ld_ref_dir = smr_eqtl_result[2]
    if len(os.listdir(_ld_ref_dir)) == 0 or len(os.listdir(_subset_vcf_dir)) == 0:
        logging.warning('Dependant files not found, did you run eqtl_data_processor?')
        return None

    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _eqtl_sample_size = glob_processor.global_config['input']['eqtl'].get('sample_size', 948)
    _genecode_file = glob_processor.global_config['input']['genecode']
    smr = rs.Smr()

    return smr.run(_working_dir,
                   gdp.Processor.VAR_ID_COL_NAME,
                   _genecode_file,
                   glob_processor.gwas_filter_file,
                   glob_processor.gwas_output_dir,
                   glob_processor.gwas_col_dict,
                   _gwas_sample_size,
                   glob_processor.eqtl_output_report,
                   glob_processor.eqtl_output_dir,
                   glob_processor.eqtl_col_dict,
                   glob_processor.config_holder.eqtl_p_threshold,
                   _subset_vcf_dir,
                   smr_eqtl_result[1],
                   _ld_ref_dir,
                   smr_eqtl_result[3],
                   tools_config_file=glob_processor.tools_config_file)


def __preprocess_and_run_ecaviar(glob_processor):
    _working_dir = os.path.join(glob_processor.tool_parent_dir, run_e.ECaviar.ECAVIAR_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _gwas_sample_size = glob_processor.global_config['input']['gwas']['sample_size']
    _eqtl_sample_size = glob_processor.global_config['input']['eqtl'].get('sample_size', 948)

    ecaviar_data_processor = edp.ECaviarDataProcessor()
    preproc_rst_dir = ecaviar_data_processor.prepare(working_dir=_working_dir,
                                                     gwas_cluster_dir=glob_processor.gwas_cluster_output_dir,
                                                     gwas_cluster_summary=glob_processor.gwas_cluster_summary,
                                                     eqtl_group_dir=glob_processor.eqtl_output_dir,
                                                     eqtl_report=glob_processor.eqtl_output_report,
                                                     ref_vcf_dir=glob_processor.ref_vcf_dir,
                                                     gwas_col_dict=glob_processor.gwas_col_dict,
                                                     eqtl_col_dict=glob_processor.eqtl_col_dict,
                                                     population=glob_processor.global_config.get('population',
                                                                                                 'EUR').upper(),
                                                     gwas_sample_size=_gwas_sample_size,
                                                     eqtl_sample_size=_eqtl_sample_size,
                                                     var_id_col_name=gdp.Processor.VAR_ID_COL_NAME)

    ecaviar = run_e.ECaviar()
    return asyncio.run(ecaviar.run(working_dir=_working_dir,
                                   candidate_data_dir=preproc_rst_dir,
                                   parallel=glob_processor.config_holder.parallel,
                                   tools_config=glob_processor.tools_config_file))


def __preprocess_and_run_twas(glob_processor):
    _working_dir = os.path.join(glob_processor.tool_parent_dir, rt.TWAS.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    if len(os.listdir(glob_processor.gwas_output_dir)) == 0:
        raise ValueError(f'Dependant GWAS files not found, did you run global_data_process?')
    _weight_pos_file = util.get_twas_ref_files(glob_processor.global_config)
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
