import argparse
import asyncio
import logging
import os
import re
import shutil
import subprocess
import sys
import warnings
from itertools import combinations
from pathlib import Path

import pandas as pd
import yaml

sys.path.append(os.path.abspath(os.path.dirname(Path(__file__).resolve())))
import constants as const

shell_command_write_vcf_header = 'tabix -H {} > {}'
shell_command_get_vcf_data_by_position = 'tabix {} chr{}:{}-{}'


def split_file_by_col_name_sqtl(working_dir, src_file_path, chrom_col, gene_col, pheno_col, readonly_cols=None, dtype=None,
                           sep='\t'):
    """

    output file names: {working_dir}/{chrom}/{gene}/{pheno}.tsv.gz
    """
    if not os.path.exists(src_file_path) or os.path.getsize(src_file_path) <= 0:
        warnings.warn(
            f'Input file {src_file_path} does not exist or size is 0, nothing to split')
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    with pd.read_table(src_file_path, sep=sep, iterator=True, chunksize=100000, header=0,
                       usecols=readonly_cols, dtype=dtype) as reader:
        for chunk in reader:
            grouped = chunk.groupby(pheno_col)
            logging.info(f"grouped: {grouped}")
            for name, group in grouped:
                chrom = group[chrom_col].iloc[0]
                gene_id = group[gene_col].iloc[0]
                subdir = os.path.join(working_dir, f'{chrom}',f'{gene_id}')
                Path(subdir).mkdir(parents=True, exist_ok=True)
                csv_file = os.path.join(subdir, f'{name}.tsv.gz')
                logging.info(f"csv_file: {csv_file}")
                # logging.info(f"*****csv_file*****: {csv_file}")
                # logging.info(f"*****name*****: {name}")
                # csv_file = os.path.join(subdir, f'{name}.tsv.gz')
                # check csv_file existence and determine write mode
                if os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
                    mode = 'a'
                    header = False
                else:
                    mode = 'w'
                    header = True
                group.to_csv(csv_file, mode=mode, sep=const.output_spliter, header=header, index=False)


def split_file_by_col_name_qtl(working_dir, src_file_path, chrom_col, pheno_col, readonly_cols=None, dtype=None,
                           sep='\t'):
    """

    output file names: {working_dir}/{chrom}/{pheno/gene}.tsv.gz
    """
    if not os.path.exists(src_file_path) or os.path.getsize(src_file_path) <= 0:
        warnings.warn(
            f'Input file {src_file_path} does not exist or size is 0, nothing to split')
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    with pd.read_table(src_file_path, sep=sep, iterator=True, chunksize=100000, header=0,
                       usecols=readonly_cols, dtype=dtype) as reader:
        for chunk in reader:
            grouped = chunk.groupby(pheno_col)
            for name, group in grouped:
                chrom = group[chrom_col].iloc[0]
                subdir = os.path.join(working_dir, f'{chrom}')
                Path(subdir).mkdir(parents=True, exist_ok=True)
                csv_file = os.path.join(subdir, f'{name}.tsv.gz')
                # logging.info(f"*****csv_file*****: {csv_file}")
                # logging.info(f"*****name*****: {name}")
                # csv_file = os.path.join(subdir, f'{name}.tsv.gz')
                # check csv_file existence and determine write mode
                if os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
                    mode = 'a'
                    header = False
                else:
                    mode = 'w'
                    header = True
                group.to_csv(csv_file, mode=mode, sep=const.output_spliter, header=header, index=False)


def filter_file_by_p_value(working_dir, input_file_path, p_value_threshold, p_value_col_name, out_file_path=None):
    result_df = None
    if not os.path.exists(input_file_path) or os.path.getsize(input_file_path) <= 0:
        warnings.warn(
            f'Input file {input_file_path} does not exist or size is 0, nothing to filter')
        return result_df
    with pd.read_table(input_file_path, sep=const.column_spliter, iterator=True, chunksize=100000, header=0) as reader:
        for chunk in reader:
            if result_df is None:
                result_df = pd.DataFrame(columns=chunk.columns)
            filter_data_frame_by_p_value(chunk, p_value_threshold, p_value_col_name, True)
            # append to result_df if there is data after filtering
            if len(chunk) > 0:
                # note that append returns a new DataFrame
                result_df = result_df.append(chunk)
    if out_file_path:
        # Write output to working_dir
        result_df.to_csv(os.path.join(working_dir, f'pval_filtered_{os.path.basename(input_file_path)}'),
                         sep=const.output_spliter, header=True, index=False)
    return result_df


def filter_data_frame_by_p_value(input_df, p_value_threshold, p_value_col_name, inplace=True):
    if not isinstance(input_df, pd.DataFrame):
        warnings.warn('input_df is not an instance of DataFrame, unable to filter')
        return None
    # drop rows with p-value > p_value_threshold
    return input_df.drop(
        input_df[(input_df[p_value_col_name] > p_value_threshold) | (input_df[p_value_col_name] == 0.0)].index,
        inplace=inplace)


def find_col_name_by_regex(regex, col_name_array):
    for col in col_name_array:
        if re.search(regex, col, re.I):
            return col
    return None


def get_file_name(file_path):
    # return file name without extension
    return file_path.split(os.sep)[-1].split('.')[0]

def get_gene_name(file_path):
    # return file name without extension
    return '.'.join(file_path.split(os.sep)[-1].split('.')[-3:-1])

def get_eqtl_gene_name(file_path):
    # print(f"get_eqtl_gene_name {file_path.split(os.sep)[-1].split('.')}")
    # return file name without extension
    if len(file_path.split(os.sep)[-1].split('.')) == 3:
        return file_path.split(os.sep)[-1].split('.')[0]
    else:
        return '.'.join(file_path.split(os.sep)[-1].split('.')[:-2])


def get_pheno_name(file_path):
    # return file name without extension
    return '.'.join(file_path.split(os.sep)[-1].split('.')[:-2])
    
def file_exists(file_path):
    if file_path is None:
        return False
    return Path(file_path).exists()


def delete_file_if_exists(file_path):
    if file_exists(file_path):
        try:
            os.remove(file_path)
            # pass
        except OSError as e:
            logging.warning(f'Fail to delete {file_path}: {e}')


def delete_dir(file_path):
    os.system(f'rm -rf {file_path}')
    # pass
    # too slow, replace by rm -rf
    # if file_exists(file_path):
    #     if Path(file_path).is_file():
    #         delete_file_if_exists(file_path)
    #     else:
    #         shutil.rmtree(file_path)


def clean_data(input_df, dup_consider_subset=None, na_consider_subset=None, keep_dup=False):
    """
    Drop na and duplicated records on the input_df

    input_df
        The input dataframe to be cleaning
    na_consider_subset
        The column label or sequence of labels used to remove na records
    dup_consider_subset
        The column label or sequence of labels used to remove duplicated records
    """
    
    
    if not isinstance(input_df, pd.DataFrame):
        warnings.warn('input_df is not an instance of DataFrame, unable to filter')
        return
    input_df.dropna(inplace=True, subset=na_consider_subset)
    # noinspection PyTypeChecker
    input_df.drop_duplicates(subset=dup_consider_subset, keep=keep_dup, inplace=True)


def clean_chunk_data(input_df, dup_consider_col, complete_dup_set_df, na_consider_subset=None):
    """
    Similar to clean_data but useful when perform cleaning on iterator chunk data
    """
    if not isinstance(input_df, pd.DataFrame):
        warnings.warn('input_df is not an instance of DataFrame, unable to filter')
        return
    input_df.dropna(inplace=True, subset=na_consider_subset)
    chunk_dup_bool_series = input_df[dup_consider_col].isin(complete_dup_set_df[dup_consider_col])
    input_df.drop(labels=input_df[chunk_dup_bool_series].index, inplace=True)


def extract_vcf_data(chromosome, target_pos_rsid_df, ref_vcf_file_path, output_vcf_dir, output_vcf_file_name,
                     target_position_col_name, target_snp_col_name, extract_step_size=500000):
    if target_pos_rsid_df.shape[0] == 0:
        logging.info(f'**** end process, target df is empty, nothing to extract')
        return
    if is_vcf_chrom_code_contains_chr(ref_vcf_file_path):
        chrom_code_in_vcf = chromosome if chromosome.startswith('chr') else f'chr{chromosome}'
    else:
        chrom_code_in_vcf = chromosome.strip('chr') if chromosome.startswith('chr') else chromosome
    # delete result file if they are already there
    output_file_full_path = f'{output_vcf_dir}/{output_vcf_file_name}'
    tmp_vcf_file_full_path = f'{output_file_full_path}.tmp'
    delete_file_if_exists(output_file_full_path)
    delete_file_if_exists(tmp_vcf_file_full_path)
    matching_snp_dir = f'{output_vcf_dir}/matching'
    report_file_name = output_vcf_file_name.replace('.vcf', '.tsv')
    output_matching_snp_report = f'{matching_snp_dir}/{report_file_name}'
    delete_file_if_exists(output_matching_snp_report)
    missing_snp_dir = f'{output_vcf_dir}/missing'
    output_missing_snp_report = f'{missing_snp_dir}/{report_file_name}'
    delete_file_if_exists(output_missing_snp_report)
    # writing header
    os.system(shell_command_write_vcf_header.format(ref_vcf_file_path, tmp_vcf_file_full_path))
    Path(matching_snp_dir).mkdir(parents=True, exist_ok=True)
    Path(missing_snp_dir).mkdir(parents=True, exist_ok=True)
    matching_positions = []
    max_pos = target_pos_rsid_df[target_position_col_name].max()
    min_pos = target_pos_rsid_df[target_position_col_name].min()
    start_pos = min_pos
    with open(tmp_vcf_file_full_path, 'a') as output_file_obj:
        while start_pos < max_pos + 1:
            # start and stop are both inclusive in tabix
            end_pos = min(start_pos + extract_step_size, max_pos)
            process = subprocess.Popen(['tabix', ref_vcf_file_path, f'{chrom_code_in_vcf}:{start_pos}-{end_pos}'],
                                       stdout=subprocess.PIPE)
            for line in iter(process.stdout.readline, b''):
                new_line = line.decode('UTF-8')
                snp_pos = int(new_line.split('\t')[1])
                if snp_pos not in target_pos_rsid_df[target_position_col_name].values:
                    continue
                matching_positions.append(snp_pos)
                split_row = new_line.split('\t')
                split_row[2] = target_pos_rsid_df[target_pos_rsid_df[target_position_col_name] == snp_pos].iloc[0].loc[
                    target_snp_col_name]
                output_file_obj.write('\t'.join(split_row))
            start_pos = start_pos + extract_step_size + 1
    # keep output vcf and matching_report_df exactly matching
    duplicate_exists = False
    pos_df = pd.DataFrame(matching_positions)
    pos_df.drop_duplicates(keep=False, inplace=True)
    # logging.info('*********pos_df********')
    # logging.info(pos_df)
    if pos_df.empty:
        logging.info(f'**** end process, pos_df is empty')
    else:
        unique_positions = pos_df[0].to_list()
        if len(matching_positions) != len(unique_positions):
            matching_positions = unique_positions
            duplicate_exists = True
        if duplicate_exists:
            if not os.path.exists(tmp_vcf_file_full_path):
                raise FileNotFoundError(f"Temporary VCF file not found: {tmp_vcf_file_full_path}")
            else:
                with open(tmp_vcf_file_full_path) as tmp_vcf, open(output_file_full_path, mode='w') as result_vcf:
                    for line in tmp_vcf:
                        if line.startswith('#'):
                            result_vcf.write(line)
                            continue
                        snp_pos = int(line.split('\t')[1])
                        if snp_pos not in matching_positions:
                            # record at snp_pos is a duplicate record
                            continue
                        result_vcf.write(line)
                os.remove(tmp_vcf_file_full_path)
        else:
            shutil.move(tmp_vcf_file_full_path, output_file_full_path)
        matching_report_df = target_pos_rsid_df[
            target_pos_rsid_df[target_position_col_name].isin(matching_positions)]
        matching_report_df.to_csv(output_matching_snp_report, sep=const.output_spliter, index=False)
        missing_df = target_pos_rsid_df[~target_pos_rsid_df[target_position_col_name].isin(matching_positions)]
        missing_df.to_csv(output_missing_snp_report, sep=const.output_spliter, index=False)
        logging.debug(f'matching {len(matching_report_df)} snps, missing {len(missing_df)} snps')
        logging.info(f'**** end process, output file is {output_file_full_path}')


def get_chromosome_number(chrom):
    if type(chrom) == str:
        return re.sub('chr', '', chrom, count=1, flags=re.IGNORECASE)
    else:
        return str(chrom)


def parse_parameters():
    parser = argparse.ArgumentParser()
    # parser.add_argument('--tools', dest='tools_list', default=['all'], type=str, nargs='*', help='coloc tools')
    # config.yml directory
    parser.add_argument('--config', dest='config_file', help='Text file with gwas, eqtl message')
    # log file path
    parser.add_argument('--log', dest='log_file', help='Log file path')
    parser.add_argument('--disable_parallel', dest='parallel', help='disable parallel', action="store_false")
    # customized parameters for each tool
    parser.add_argument('--tools_config', dest='tools_config', help='customized parameters for each tool')
    parser.add_argument('--no_report', dest='no_report', help='turn off generating report', action="store_true")
    # parser.add_argument('--target_loci', dest='target_loci', help='Target GWAS loci range') # chr1:xxxx-xxxx
    parse_args = parser.parse_args()
    return parse_args


def read_config(config_file):
    logging.info(f'read config file: {config_file}')
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    config_dir = Path(config_file).resolve().parent.parent
    
    def update_path(value):
        if isinstance(value, str):
            if ('/' in value or '\\' in value) and not value.startswith(("\t", "\n", "\r")):
                return str(config_dir / value) if not Path(value).is_absolute() else value
        return value
    def recursive_update(d):
        if isinstance(d, dict):
            return {key: recursive_update(update_path(value)) for key, value in d.items()}
        elif isinstance(d, list):
            return [recursive_update(update_path(item)) for item in d]
        else:
            return d
    config = recursive_update(config)
    return config


def drop_diff_chrom_rows(ref_snp, snp_col_name, target_df, chrom_col_name):
    
    
    chrom = target_df[target_df[snp_col_name] == ref_snp].iloc[0][chrom_col_name]
    diff_chr_series = target_df[chrom_col_name] != chrom
    target_df.drop(labels=target_df[diff_chr_series].index, inplace=True)


def drop_non_intersect_rows(df1, df1_ref_col_name, df2, df2_ref_col_name):

    df1_intersect_series = df1[df1_ref_col_name].isin(df2[df2_ref_col_name])
    df1.drop(labels=df1[~df1_intersect_series].index, inplace=True)

    df2_intersect_series = df2[df2_ref_col_name].isin(df1[df1_ref_col_name])
    df2.drop(labels=df2[~df2_intersect_series].index, inplace=True)


def adjust_allele_order(gwas_df, ea_col_name, oa_col_name, gwas_chrom_col_name, gwas_pos_col_name,
                        ref_df, ref_df_chrom_col_name='chromosome', ref_df_pos_col_name='position',
                        ref_df_alt_allele_col_name='alt', ref_df_ref_allele_col_name='ref',
                        gbeta_col_name=None, geaf_col_name=None, gz_col_name=None,
                        drop_src_df_non_intersect_items=True,
                        drop_ref_df_non_intersect_items=True):
    logging.info(f"adjust_allele_order")
    """
        Adjust gwas_df beta/zscore sign and eaf value according a reference df allele order,
        the reference df can be an eqtl df or vcf df

        returns None

        gwas_df
            the ref gwas df
        ea_col_name
            the effect_allele column name in gwas_df
        oa_col_name
            the non_effect_allele column name in gwas_df
        gwas_chrom_col_name
            the chromosome column name in gwas_df
        gwas_pos_col_name
            the position column name in gwas_df
        ref_df
            the target reference df, can be eqtl df or vcf df
        ref_df_chrom_col_name
            the chromosome column name in ref_df
        ref_df_pos_col_name
            the position column name in target ref_df
        ref_df_alt_allele_col_name
            the alt column name in target ref_df
        ref_df_ref_allele_col_name
            the ref column name in target ref_df
        gbeta_col_name
            the beta column name in target gwas_df, sign will be adjusted according eqtl_df allele order
        geaf_col_name
            the eaf column name in target gwas_df, value will be adjusted according eqtl_df allele order
        gz_col_name
            the zscore column name in target gwas_df, sign will be adjusted according eqtl_df allele order
        drop_ref_df_non_intersect_items
            whether to drop non-intersect records in eqtl_df
        """
    merged_alt_col_name = f'{ref_df_alt_allele_col_name}_qtl' \
        if ref_df_alt_allele_col_name in [ea_col_name, oa_col_name] else ref_df_alt_allele_col_name
    merged_ref_col_name = f'{ref_df_ref_allele_col_name}_qtl' \
        if ref_df_ref_allele_col_name in [ea_col_name, oa_col_name] else ref_df_ref_allele_col_name
    merged = pd.merge(left=gwas_df[[gwas_chrom_col_name, gwas_pos_col_name, ea_col_name, oa_col_name, ]],
                      right=ref_df[[ref_df_chrom_col_name, ref_df_pos_col_name, ref_df_alt_allele_col_name, ref_df_ref_allele_col_name]],
                      left_on=[gwas_chrom_col_name, gwas_pos_col_name],
                      right_on=[ref_df_chrom_col_name, ref_df_pos_col_name], 
                      how='left', suffixes=(None, '_qtl'))
    print(f"merged.shape: {merged.shape}")
    print(merged)
    print(f"gwas_df.shape: {gwas_df.shape}")
    print(gwas_df)
    # logging.info(f"merged_alt_col_name: {merged_alt_col_name}")
    # logging.info(f"merged_ref_col_name: {merged_ref_col_name}")
    # logging.info(f"len(merged): {len(merged)}")
    # logging.info(f"merged: {merged}")
    # logging.info(f"ea_col_name: {ea_col_name}")
    # logging.info(f"oa_col_name: {oa_col_name}")

    # align index to gwas_df index
    merged.index = gwas_df.index.copy()
    merged.drop(columns=[col for col in merged.columns if
                         col not in [ea_col_name, oa_col_name, merged_alt_col_name, merged_ref_col_name]], inplace=True)
    # set element ignores order, so:
    # (gwas.ea = eqtl.alt & gwas.oa = eqtl.ref) or (gwas.ea = eqtl.ref & gwas.oa = eqtl.alt)
    # will result true in eq_bool_series
    # set is unhashable, can not be used as CategoricalDtype, use frozenset instead.
    # category dtype will drastically reduce memory use.


    # 确保列为字符串类型，防止 Categorical 类型引发错误
    merged[ea_col_name] = merged[ea_col_name].astype(str)
    merged[oa_col_name] = merged[oa_col_name].astype(str)
    merged[merged_alt_col_name] = merged[merged_alt_col_name].astype(str)
    merged[merged_ref_col_name] = merged[merged_ref_col_name].astype(str)

    # allele_cat = pd.CategoricalDtype([frozenset(e) for e in combinations(const.SNP_ALLELE, 2)])
    gwas_alleles = pd.Series([frozenset(e) for e in zip(merged[ea_col_name], merged[oa_col_name])],
                             index=gwas_df.index.copy(), 
                            #  dtype=allele_cat
                             )
    ref_panel_alleles = pd.Series([frozenset(e) for e in zip(merged[merged_alt_col_name], merged[merged_ref_col_name])],
                                  index=gwas_df.index.copy(), 
                                #   dtype=allele_cat
                                  )
    
    eq_bool_series = gwas_alleles == ref_panel_alleles
    del gwas_alleles, ref_panel_alleles
    gwas_flipped_bool_series = merged[ea_col_name] != merged[merged_alt_col_name]
    eq_flipped_list = (eq_bool_series & gwas_flipped_bool_series).to_list()
    del gwas_flipped_bool_series
    logging.info(f"gwas_df[ea_col_name] dtype: {gwas_df[ea_col_name].dtype}")
    logging.info(f"gwas_df[oa_col_name] dtype: {gwas_df[oa_col_name].dtype}")
    logging.info(f"merged[merged_alt_col_name] dtype: {merged[merged_alt_col_name].dtype}")
    logging.info(f"merged[merged_ref_col_name] dtype: {merged[merged_ref_col_name].dtype}")

    gwas_df.loc[eq_flipped_list, ea_col_name] = merged.loc[eq_flipped_list, merged_alt_col_name]
    gwas_df.loc[eq_flipped_list, oa_col_name] = merged.loc[eq_flipped_list, merged_ref_col_name]
    del merged
    if gbeta_col_name is not None:
        gwas_df.loc[eq_flipped_list, gbeta_col_name] = - gwas_df.loc[eq_flipped_list, gbeta_col_name]
    if gz_col_name is not None:
        gwas_df.loc[eq_flipped_list, gz_col_name] = - gwas_df.loc[eq_flipped_list, gz_col_name]
    if geaf_col_name is not None:
        gwas_df.loc[eq_flipped_list, geaf_col_name] = 1 - gwas_df.loc[eq_flipped_list, geaf_col_name]
    if drop_src_df_non_intersect_items:
        gwas_df.drop(index=gwas_df[~eq_bool_series].index, inplace=True)
    if drop_ref_df_non_intersect_items:
        ref_df_inter_bool_series = ref_df[ref_df_chrom_col_name].isin(gwas_df[gwas_chrom_col_name]) & ref_df[
            ref_df_pos_col_name].isin(gwas_df[gwas_pos_col_name])
        ref_df.drop(labels=ref_df[~ref_df_inter_bool_series].index, inplace=True)
    logging.info(f"Adjust_allele_order Done")


def drop_indel_snp(target_df, allele_col_name1, allele_col_name2):
    """
        Drop INDEL SNPs from target_df inplace

        returns None

        target_df
            the target df
        allele_col_name1
            the effect_allele/alt or  other_allele/ref column name in target_df
        allele_col_name2
            the other_allele/ref or effect_allele/alt column name in target_df
        """
    indel_bool_series = (target_df[allele_col_name1].str.len() != 1) | (target_df[allele_col_name2].str.len() != 1)
    target_df.drop(labels=target_df[indel_bool_series].index, inplace=True)


def check_path_exist_and_has_size(file, raise_error=True):
    if not os.path.exists(file) or os.path.getsize(file) <= 0:
        if raise_error:
            raise ValueError(f'file: {file} ,Dependant files not found')
        else:
            logging.warning(f'{file} no such file or directory')


def check_file_or_path_exist(file_path, raise_error=True):
    if not file_exists(file_path):
        if raise_error:
            raise ValueError(f'{file_path} no such file or directory')
        else:
            logging.warning(f'{file_path} no such file or directory')


def mapping_var_id_to_rsid(result_df, result_df_var_id_col_name,
                           result_df_gene_id_col_name, gwas_preprocessed_file=None,
                           ref_var_id_col_name=None, gwas_col_dict=None, qtl_col_dict=None):
    if result_df is None or result_df.shape[0] == 0:
        return result_df
    if gwas_preprocessed_file is not None and gwas_col_dict.get('snp') is not None and ref_var_id_col_name is not None:
        mapping_df = pd.read_table(gwas_preprocessed_file, usecols=[ref_var_id_col_name, gwas_col_dict['snp']])
        mapping_df.rename({gwas_col_dict['snp']: 'rsid'}, axis='columns', inplace=True)
        # Now 2 columns left in mapping_df: ref_var_id_col_name, rsid
        result_df = pd.merge(left=result_df, right=mapping_df,
                             left_on=result_df_var_id_col_name, right_on=ref_var_id_col_name,
                             how='left')
        # Drop ref_var_id_col_name from result_df if we introduce extra ref_var_id_col_name
        if result_df_var_id_col_name != ref_var_id_col_name:
            result_df.drop(columns=ref_var_id_col_name, inplace=True)
    elif qtl_col_dict is not None and qtl_col_dict.get('snp') is not None and ref_var_id_col_name is not None:
        merged_dfs = []
        merged_genes = set()
        for indx, row in result_df.iterrows():
            gene_id = row.loc[result_df_gene_id_col_name]
            eqtl_gene_file = row.loc['eqtl_path']
            if (gene_id in merged_genes) or (not file_exists(eqtl_gene_file)):
                continue

            mapping_df = pd.read_table(eqtl_gene_file, usecols=[ref_var_id_col_name, qtl_col_dict['snp']])
            mapping_df.rename({qtl_col_dict['snp']: 'rsid'}, axis='columns', inplace=True)
            # merged has 2 or 3 columns:
            # result_df_var_id_col_name (if not equal to ref_var_id_col_name), ref_var_id_col_name, rsid
            merged = pd.merge(left=result_df[result_df_var_id_col_name], right=mapping_df,
                              left_on=result_df_var_id_col_name, right_on=ref_var_id_col_name, how='inner')
            if result_df_var_id_col_name != ref_var_id_col_name:
                merged.drop(columns=ref_var_id_col_name, inplace=True)
            merged_dfs.append(merged)
            merged_genes.add(gene_id)
        del merged_genes
        if len(merged_dfs) == 0:
            return result_df
        merged_df = pd.concat(merged_dfs)
        del merged_dfs
        # merged_df has 2 columns: result_df_var_id_col_name, rsid
        # merged_df may have duplicates since diff eqtl gene file may contain the same subset of snp
        merged_df.drop_duplicates(inplace=True)
        result_df = pd.merge(left=result_df, right=merged_df,
                             left_on=result_df_var_id_col_name, right_on=result_df_var_id_col_name,
                             how='left')
    return result_df


def get_predixcan_ref_files(global_config):
    prediction_source_dir = global_config['input']['prediction_dir']
    eqtl_issue_name = global_config['input']['qtl']['biological_context']
    model_db_path = None
    prediction_snp_covariance_path = None
    for filename in os.listdir(prediction_source_dir):
        if model_db_path is not None and prediction_snp_covariance_path is not None:
            break
        if re.search(rf'^(.*){eqtl_issue_name}(.*)\.db$', filename):
            model_db_path = os.path.join(prediction_source_dir, filename)
        elif re.search(rf'^(.*){eqtl_issue_name}(.*)\.txt.gz$', filename):
            prediction_snp_covariance_path = os.path.join(prediction_source_dir, filename)
    return model_db_path, prediction_snp_covariance_path


def get_twas_ref_files(global_config, filtered_pos=False):
    twas_model_dir = global_config['input']['twas_model_dir']
    eqtl_issue_name = global_config['input']['qtl']['biological_context']
    pos_name = None
    no_filter_pos_name = None
    for filename in os.listdir(twas_model_dir):
        if pos_name is not None and no_filter_pos_name is not None:
            break
        if re.search(rf'^(.*){eqtl_issue_name}(.*)\.pos$', filename):
            if 'nofilter' in filename:
                no_filter_pos_name = filename
            else:
                pos_name = filename
    logging.info(f"pos_name: {pos_name}")
    logging.info(f"path: {os.path.join(twas_model_dir, pos_name)}")
    if no_filter_pos_name is not None:
        pos_file_name = no_filter_pos_name
    else:
        pos_file_name = pos_name
    # pos_file_name = pos_name if filtered_pos else no_filter_pos_name
    logging.info(f"pos_file_name: {pos_file_name}")
    logging.info(f"filtered_pos: {filtered_pos}")
    return os.path.join(twas_model_dir, pos_file_name)


def is_vcf_chrom_code_contains_chr(vcf_file):
    if not (Path(f'{vcf_file}.tbi').exists() or Path(f'{vcf_file}.csi').exists()):
        os.system(f'tabix -f -p vcf {vcf_file}')
    process = subprocess.Popen(['tabix', '-l', vcf_file],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        return False
    return 'chr' in stdout.decode('UTF-8')


async def async_run_cmd(cmd, logging_stdout=False):
    
    
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=sys.stdout)

    stdout, stderr = await proc.communicate()

    logging.info(f'[{cmd!r} exited with {proc.returncode}]')
    if stdout and logging_stdout:
        logging.info(f'[stdout]\n{stdout.decode()}')
    if stderr:
        logging.info(f'[stderr]\n{stderr.decode()}')


async def gather_with_limit(limit, *coros):
    """
    asyncio.gather alternative with limit number of concurrent coroutines, mainly used to avoid OOM

    limit
        Number of concurrent coroutines
    coros
        coroutines, NOT tasks
    """
    
    
    semaphore = asyncio.Semaphore(limit)

    async def sem_coro(coro):
        async with semaphore:
            return await coro

    return await asyncio.gather(*(sem_coro(c) for c in coros), return_exceptions=True)


def run_logging_command(command):
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=sys.stdout)

    while True:
        line = p.stdout.readline().decode("utf-8")
        if not line:
            break
        logging.info(line)


def cleanup_output(tools_output_base):
    if file_exists(tools_output_base):
        for tool_output_dir in os.listdir(tools_output_base):
            if tool_output_dir == 'rank':
                continue
            tool_output_dir_full_path = f'{tools_output_base}/{tool_output_dir}'
            if Path(tool_output_dir_full_path).is_dir():
                for tool_output_sub_dir in os.listdir(tool_output_dir_full_path):
                    if tool_output_sub_dir == 'analyzed':
                        continue
                    tool_output_sub_out_full_path = os.path.join(tool_output_dir_full_path, tool_output_sub_dir)
                    if Path(tool_output_sub_out_full_path).is_dir():
                        delete_dir(tool_output_sub_out_full_path)
                    else:
                        delete_file_if_exists(tool_output_sub_out_full_path)
            else:
                delete_file_if_exists(tool_output_dir_full_path)


# remove nan rows and columns in ld_file, return removed column names
def remove_nan_from_ld(ld_file, header):
    # remove last empty column generated by plink
    ld_df = pd.read_csv(ld_file, sep=' ', header=None)
    ld_row_size = len(ld_df)
    ld_col_size = len(ld_df.columns)
    if (ld_col_size - ld_row_size == 1) and (ld_df[len(ld_df.columns) - 1].isnull().sum() == len(ld_df)):
        ld_df.drop(columns=ld_col_size - 1, inplace=True)
    ld_df.columns = header

    nan_df = ld_df.columns[ld_df.isna().sum() == len(ld_df)]
    nan_cols = nan_df.tolist()
    # remove nan columns
    ld_df.drop(nan_cols, axis=1, inplace=True)
    # remove nan rows
    ld_df.dropna(how='all', axis=0, inplace=True)
    ld_df.to_csv(ld_file, mode='w', sep=' ', index=False, header=False)
    return nan_cols


def get_tools_params(tool_name, tools_config_file, param_prefix='--', params_without_value=[]):
    params_str = ''
    param_dict = get_tools_params_dict(tool_name, tools_config_file)
    logging.info(f"param_dict: {param_dict}")
    if param_dict is not None:
        for pk, pv in param_dict.items():
            if pv is None:
                if pk in params_without_value:
                    params_str += f'{param_prefix}{pk} '
            else:
                params_str += f'{param_prefix}{pk} {pv} '
    else:
        logging.info(f'No parameter config for {tool_name}, use default parameters')
    return params_str


def get_tools_params_dict(tool_name, tool_config_file):
    params = None
    if tool_config_file and file_exists(tool_config_file):
        with open(tool_config_file, 'r') as cfg_file:
            params = yaml.safe_load(cfg_file)
            if params is not None and tool_name in params:
                params = params[tool_name]
            else:
                logging.info(f'No parameter config for {tool_name}, use default parameters')
                params = None
    else:
        logging.info(f'No parameter config for {tool_name}, use default parameters')
    return {} if params is None else params