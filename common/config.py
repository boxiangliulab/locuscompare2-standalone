from pathlib import Path
import os
import common.constants
from common import coloc_utils as utils


class ConfigHolder:

    def __init__(self, single_config_file=common.constants.default_config, study=common.constants.default_study,
                 parallel=False):
        if study is None:
            study = common.constants.default_study
        if single_config_file is None:
            single_config_file = common.constants.default_config

        self.global_config = utils.read_config(single_config_file)
        self.root_work_dir = self.global_config['working_dir']
        self.output_preprocessed_dir = os.path.join(self.root_work_dir, 'preprocessed')
        self.output_processed_dir = os.path.join(self.root_work_dir, 'processed')
        Path(self.output_preprocessed_dir).mkdir(exist_ok=True, parents=True)
        Path(self.output_processed_dir).mkdir(exist_ok=True, parents=True)

        self.gwas_trait = self.global_config['input']['gwas']['trait']
        self.eqtl_tissue = self.global_config['input']['eqtl']['tissue']
        self.population = self.global_config['population']

        self.study_dir = os.path.join(self.output_processed_dir, study)

        self.tool_parent_dir = os.path.join(self.study_dir, self.gwas_trait, self.eqtl_tissue, self.population)
        Path(self.tool_parent_dir).mkdir(exist_ok=True, parents=True)

        self.report_path = os.path.join(self.study_dir)

        self.gwas_col_dict = self.global_config['input']['gwas']['col_name_mapping']
        self.eqtl_col_dict = self.global_config['input']['eqtl']['col_name_mapping']

        # eqtl output path
        self.eqtl_output_dir = os.path.join(self.output_preprocessed_dir, 'eqtl', self.eqtl_tissue, 'grouped')
        self.eqtl_output_report = os.path.join(self.output_preprocessed_dir, 'eqtl', self.eqtl_tissue,
                                               'filtered_gene.tsv.gz')

        # gwas output path
        self.gwas_preprocessed_dir = os.path.join(self.output_preprocessed_dir, 'gwas', self.gwas_trait)
        self.gwas_preprocessed_file = os.path.join(self.gwas_preprocessed_dir, 'preprocessed.tsv.gz')
        self.gwas_output_dir = os.path.join(self.output_preprocessed_dir, 'gwas', self.gwas_trait, 'grouped')
        self.gwas_filter_file = os.path.join(self.gwas_preprocessed_dir, 'pval_filtered.tsv.gz')
        self.gwas_cluster_output_dir = os.path.join(self.gwas_preprocessed_dir, 'clustered')
        self.gwas_cluster_summary = os.path.join(self.gwas_preprocessed_dir, 'cluster_summary.tsv.gz')

        # vcf input output path
        self.vcf_output_dir = os.path.join(self.gwas_preprocessed_dir, 'vcf')
        self.ref_vcf_dir = self.global_config['input']['vcf']

        self.parallel = parallel

        # p-val threshold
        self.gwas_p_threshold = self.global_config['p-value_threshold']['gwas']
        if self.gwas_p_threshold <= 0:
            self.gwas_p_threshold = 5.0E-8
        elif self.gwas_p_threshold > 1.0E-7:
            self.gwas_p_threshold = 1.0E-7
        self.eqtl_p_threshold = self.global_config['p-value_threshold']['eqtl']
        if self.eqtl_p_threshold <= 0:
            self.eqtl_p_threshold = 1.0E-6
        elif self.eqtl_p_threshold > 1.0E-5:
            self.eqtl_p_threshold = 1.0E-5
