import os
import sys
from pathlib import Path
import pandas as pd
import subprocess
import argparse
from datetime import datetime


# tools -> dap torus perl bgzip gzip qtltools
#   genotype file the third columns format decide the third columns format of eqtl annotations file finnally

# step 1 covariate_file no need zip
def converting_to_sbams_format(eqtl_path, phenotype_file, genotype_file, covariate_file, tissue, current_perl_path):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    print(f'converting_to_sbams_format')

    output_dir = f'{eqtl_path}/{tissue}'
    shell_command_perl_process_execute = 'perl {} -e {} -g {} -c {} -t {} -a {}'
    shell_command_bash_perl_cmd = 'bash {}'

    cmd_str_print_and_os_system(
        shell_command_perl_process_execute.format(f'{current_perl_path}/process.pl', phenotype_file, genotype_file,
                                                  covariate_file, output_dir, f'{current_perl_path}/assemble.pl'))
    cmd_str_print_and_os_system(shell_command_bash_perl_cmd.format(f'{eqtl_path}/{tissue}.assemble.cmd'))

    return output_dir


def get_chroms_from_bed(vcf_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    if not (Path(f'{vcf_file}.tbi').exists() or Path(f'{vcf_file}.csi').exists()):
        os.system(f'tabix -f -p bed {vcf_file}')
    process = subprocess.Popen(['tabix', '-l', vcf_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        return []
    return stdout.decode("UTF-8").splitlines()


# step 2 covariate_file need zip
def estimate_priors_for_finemapping(eqtl_path, phenotype_file, genotype_file, covariate_file):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    print(f'estimate_priors_for_finemapping')
    fastqtl_dir = f'{eqtl_path}/fastqtl'
    Path(fastqtl_dir).mkdir(parents=True, exist_ok=True)

    # get fastqtl file
    # foreach 22 chrom to get fastqtl file,then concat to one file
    torus_input_file_list = []
    torus_input_file = f'{eqtl_path}/torus_input.tsv'
    chroms = get_chroms_from_bed(phenotype_file)
    for chr_num in chroms:
        qtltools_file = f'{fastqtl_dir}/priors_{chr_num}.txt'
        shell_command_qtltools_execute = 'QTLtools cis --vcf {} --bed {} --cov {} --nominal 0.01 --region {} --normal --out {} --std-err'
        cmd_str_print_and_os_system(
            shell_command_qtltools_execute.format(genotype_file, phenotype_file, f'{covariate_file}.gz',
                                                  chr_num,
                                                  qtltools_file))
        # input data format
        # column 1: gene name
        # column 2: SNP name
        # column 3: distance to TSS
        # column 4: p-value
        # column 5: beta-hat
        df_qtltools_file = pd.read_csv(qtltools_file, sep=r'\s+', header=None)
        torus_input_file_list.append(df_qtltools_file[[0, 7, 6, 11, 13, 14, 9]])
    pd.concat(torus_input_file_list).to_csv(
        torus_input_file,
        sep='\t',
        index=False,
        header=False)

    shell_compress_file = 'gzip -k -f {}'
    os.system(shell_compress_file.format(torus_input_file))

    # torus to prior file
    torus_output_dir = f'{eqtl_path}/torus_output'
    shell_command_torus_execute = 'torus -d {} --fastqtl -dump_prior {}'
    cmd_str_print_and_os_system(shell_command_torus_execute.format(f'{torus_input_file}.gz', torus_output_dir))
    return torus_output_dir


# step 3
def annotations_finemapping_analysis(eqtl_path, sbams_files_dir, priors_files_dir):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    print('annotations_finemapping_analysis')
    # dap-g to finemapping
    shell_command_dap_execute = 'dap-g -d {} -p {} -ld_control 0.5 --all -t 4 > {}'
    dag_output_dir = f'{eqtl_path}/postfix'
    Path(dag_output_dir).mkdir(parents=True, exist_ok=True)
    for priors_file_name in os.listdir(priors_files_dir):
        if priors_file_name.endswith('.prior'):
            gene_name = priors_file_name.replace('.prior', '')
            sbams_file_name = f'{gene_name}.sbams.dat'
            sbams_file = f'{sbams_files_dir}/{sbams_file_name}'
            priors_file = f'{priors_files_dir}/{priors_file_name}'
            if Path(sbams_file).exists():
                dag_output_file = f'{dag_output_dir}/{gene_name}.postfix'

                cmd_str_print_and_os_system(shell_command_dap_execute.format(sbams_file,
                                                                             priors_file,
                                                                             dag_output_file))
                print(f'postfix:{dag_output_file}')
            else:
                print(f'prior file:{priors_file} has no mapping sbams file:{sbams_file}')
    return dag_output_dir


# step 4
def derive_annotations_based_on_individual_data(eqtl_path, dap_rst_dir, snp_vcf_file, tissue, current_perl_path,
                                                output_path):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    print(f'derive_annotations_based_on_individual_data')
    # perl script to annotations vcf file
    shell_command_annotation_execute = 'perl {} -dir {} -vcf {} -tissue {} | bgzip > {}'
    cmd_str_print_and_os_system(
        shell_command_annotation_execute.format(f'{current_perl_path}/summarize_dap2enloc.pl', dap_rst_dir,
                                                snp_vcf_file, tissue, output_path))
    return output_path


def prepare(parse_args):
    start_time = datetime.now()
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    print(f'prepare start at: {start_time}')

    eqtl_path = parse_args.target_dir
    genotype_file = parse_args.genotype_file
    phenotype_file = parse_args.phenotype_file
    covariate_file = parse_args.covariate_file
    tissue_name = parse_args.tissue_name
    current_perl_path = parse_args.perl_dir
    tissue = tissue_name if tissue_name else 'tissue'

    Path(eqtl_path).mkdir(parents=True, exist_ok=True)
    sbams_files_dir = converting_to_sbams_format(eqtl_path, phenotype_file, genotype_file, covariate_file,
                                                 tissue, current_perl_path)
    priors_files_dir = estimate_priors_for_finemapping(eqtl_path, phenotype_file, genotype_file,
                                                       covariate_file)
    dag_output_dir = annotations_finemapping_analysis(eqtl_path, sbams_files_dir, priors_files_dir)
    eqtl_processor_data_file = derive_annotations_based_on_individual_data(eqtl_path, dag_output_dir,
                                                                           genotype_file, tissue, current_perl_path,
                                                                           parse_args.output_path)
    print(f'eqtl annotation vcf file:{eqtl_processor_data_file}')
    print(f'prepare complete at: {datetime.now()},duration: {datetime.now() - start_time}')


def parse_parameters():
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    parser = argparse.ArgumentParser()
    parser.add_argument('--d', dest="target_dir", help="target annotation file path")
    parser.add_argument('--g', dest="genotype_file", help="Text file with genotype vcf zip file")
    parser.add_argument('--p', dest="phenotype_file", help="Text file with phentype bed zip file")
    # covariate_file require both gzip and unzip version
    parser.add_argument('--c', dest="covariate_file",
                        help="Text file with covariate txt file, but covariate zip file is also required in the current directory")
    parser.add_argument('--perl_dir', dest="perl_dir", help="The directory of the perl script")
    parser.add_argument('--t', dest="tissue_name", default='tissue', help="the tissue name")
    parser.add_argument('--o', dest="output_path", help="The output file path")
    parse_args = parser.parse_args()
    return parse_args


def cmd_str_print_and_os_system(cmd_str):
    print(cmd_str)
    os.system(cmd_str)


if __name__ == '__main__':
    args = parse_parameters()
    prepare(args)
