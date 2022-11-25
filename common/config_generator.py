import yaml
import os
from pathlib import Path
import argparse


def generate_gwas_eqtl_config(output_config_dir, config_template_path=None, gwas_config_file_path=None,
                              eqtl_config_file_path=None):
    default_resource_dir_path = f'{Path(__file__).parent.parent.resolve()}/resource'
    input_gwas_config_file_path = gwas_config_file_path if gwas_config_file_path else os.path.abspath(
        f'{default_resource_dir_path}/gwas_config.yml')
    input_eqtl_config_file_path = eqtl_config_file_path if eqtl_config_file_path else os.path.abspath(
        f'{default_resource_dir_path}/eqtl_config.yml')
    input_config_template_path = config_template_path if config_template_path else os.path.abspath(
        f'{default_resource_dir_path}/config_template.yml')

    with open(input_gwas_config_file_path, 'r') as gfile:
        gwas_configs = yaml.safe_load(gfile)

    with open(input_eqtl_config_file_path, 'r') as efile:
        eqtl_configs = yaml.safe_load(efile)

    with open(input_config_template_path, 'r') as template_file:
        config_template = yaml.safe_load(template_file)

    if not Path(output_config_dir).exists():
        os.mkdir(output_config_dir)

    for g_config in gwas_configs:
        for e_config in eqtl_configs:
            output_config_file_name = f'{g_config["trait"]}_{e_config["tissue"]}_config.yml'
            config_template['input']['gwas'] = g_config
            config_template['input']['eqtl'] = e_config
            config_yaml_text = yaml.dump(config_template)
            with open(f'{output_config_dir}/{output_config_file_name}', 'w') as output_config_file:
                output_config_file.write(config_yaml_text)


def parse_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', dest='output_dir', default='', type=str, nargs='*',
                        help="The generated config files output directory")
    parser.add_argument('--gb_temp', dest='global_template', default=None, type=str, nargs='*',
                        help="Global config template")
    parser.add_argument('--gw_temp', dest='gwas_template', default=None, type=str, nargs='*',
                        help="GWAS config template")
    parser.add_argument('--eqtl_temp', dest='eqtl_template', default=None, type=str, nargs='*',
                        help="eQTL config template")
    return parser.parse_args()


if __name__ == '__main__':
    parse_args = parse_parameters()
    generate_gwas_eqtl_config(parse_args.output_dir[0], parse_args.global_template, parse_args.gwas_template,
                              parse_args.eqtl_template)
