import yaml
import os
import sys
from pathlib import Path
import argparse
import pandas as pd


def generate_gwas_eqtl_config(output_config_dir, gwas_dir, eqtl_dir):
    
    
    default_resource_dir_path = f'{Path(__file__).parent.parent.resolve()}/resource'
    ref_resource_dir_path = f'{Path(__file__).parent.parent.resolve()}/doc'
    input_gwas_config_file_path = f'{default_resource_dir_path}/gwas_config.yml'
    input_eqtl_config_file_path = f'{default_resource_dir_path}/eqtl_config.yml'
    input_config_template_path = f'{default_resource_dir_path}/config_template.yml'
    silver_std_gwas_ref_file = f'{ref_resource_dir_path}/silver_std_traits.csv'
    eqtl_ref_file = f'{ref_resource_dir_path}/GTEx_v8_tissue_list.tsv'

    with open(input_gwas_config_file_path, 'r') as gfile:
        gwas_configs = yaml.safe_load(gfile)

    with open(input_eqtl_config_file_path, 'r') as efile:
        eqtl_configs = yaml.safe_load(efile)

    with open(input_config_template_path, 'r') as template_file:
        config_template = yaml.safe_load(template_file)

    if not Path(output_config_dir).exists():
        os.mkdir(output_config_dir)

    gwas_ref_df = pd.read_csv(silver_std_gwas_ref_file, usecols=['ACCESSION_ID', 'trait', 'N_num'])
    eqtl_ref_df = pd.read_table(eqtl_ref_file)

    for gwas_file in os.listdir(gwas_dir):
        if gwas_file.endswith('.gz'):
            accession_id = gwas_file.split('-')[1]
            gwas_record_df = gwas_ref_df[gwas_ref_df['ACCESSION_ID'] == accession_id]
            gwas_sample_size = int(gwas_record_df.iloc[0]['N_num'])
            trait_name = gwas_record_df.iloc[0]['trait']
            gwas_path = f'{gwas_dir}/{gwas_file}'
            gwas_configs[0]['file'] = gwas_path
            gwas_configs[0]['trait'] = trait_name
            gwas_configs[0]['sample_size'] = gwas_sample_size
            output_sub_dir = f'{output_config_dir}/{trait_name}'
            Path(output_sub_dir).mkdir(parents=True, exist_ok=True)

            for eqtl_file in os.listdir(eqtl_dir):
                if eqtl_file.endswith('.gz'):
                    tissue_name = eqtl_file.split('.')[0]
                    eqtl_record_df = eqtl_ref_df[eqtl_ref_df['tissue_name'] == tissue_name]
                    eqtl_sample_size = int(eqtl_record_df.iloc[0]['sample_size'])
                    eqtl_path = f'{eqtl_dir}/{eqtl_file}'

                    # generate config.yml
                    output_config_file_name = f'{trait_name}_{tissue_name}_config.yml'
                    eqtl_configs[0]['file'] = eqtl_path
                    eqtl_configs[0]['tissue'] = tissue_name
                    eqtl_configs[0]['sample_size'] = eqtl_sample_size
                    config_template['input']['gwas'] = gwas_configs[0]
                    config_template['input']['eqtl'] = eqtl_configs[0]
                    config_yml_text = yaml.dump(config_template)
                    with open(f'{output_sub_dir}/{output_config_file_name}', 'w') as output_config_file:
                        output_config_file.write(config_yml_text)


def parse_parameters():
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', dest='output_dir', default='', type=str, nargs='*',
                        help="The generated config files output directory")
    parser.add_argument('--gwas', dest='input_gwas_dir', default=None, type=str, nargs='*',
                        help="Global config template")
    parser.add_argument('--eqtl', dest='input_eqtl_dir', default=None, type=str, nargs='*',
                        help="GWAS config template")
    return parser.parse_args()


if __name__ == '__main__':
    parse_args = parse_parameters()
    generate_gwas_eqtl_config(parse_args.output_dir[0], parse_args.input_gwas_dir[0], parse_args.input_eqtl_dir[0])
