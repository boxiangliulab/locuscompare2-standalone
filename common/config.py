from pathlib import Path
import os
import sys
import common.constants
from common import utils
import logging


class ConfigHolder:

    def __init__(self, single_config_file=common.constants.default_config, study=common.constants.default_study,
                 parallel=False, tools_config_file=None):
        if study is None:
            study = common.constants.default_study
        if single_config_file is None:
            single_config_file = common.constants.default_config

        self.global_config = utils.read_config(single_config_file)
        self.root_work_dir = self.global_config['working_dir']
        print(f"global_config: {self.global_config}")
        self.qtl_type = self.global_config['input']['qtl']['qtl_type']
        self.tools = self.global_config['tools']
        self.output_preprocessed_dir = os.path.join(self.root_work_dir, 'preprocessed')
        self.output_processed_dir = os.path.join(self.root_work_dir, 'processed')
        Path(self.output_preprocessed_dir).mkdir(exist_ok=True, parents=True)
        Path(self.output_processed_dir).mkdir(exist_ok=True, parents=True)

        self.gwas_trait = self.global_config['input']['gwas']['trait'] 
        # tissue or celltype for QTLs
        self.biological_context = self.global_config['input']['qtl']['biological_context'] 
        self.population = self.global_config['population']

        self.study_dir = os.path.join(self.output_processed_dir, study)

        self.gwas_genomic_window = self.global_config['input']['gwas']['genomic_window']
        self.gwas_window_size = self.global_config['input']['gwas']['window_size']
        self.gwas_LD_r2_filter = self.global_config['input']['gwas']['LD_r2_filter']
        self.gwas_LD_additional_expansion = self.global_config['input']['gwas']['LD_additional_expansion']

        self.qtl_genomic_window = self.global_config['input']['qtl']['genomic_window']
        self.qtl_window_size = self.global_config['input']['qtl']['window_size']
        self.qtl_LD_r2_filter = self.global_config['input']['qtl']['LD_r2_filter']
        self.qtl_LD_additional_expansion = self.global_config['input']['qtl']['LD_additional_expansion']

        if self.gwas_genomic_window == self.qtl_genomic_window == 'LD_based_window':
            self.global_LD_based_window = True
        else:
            self.global_LD_based_window = False

        self.target_loci = self.global_config['input']['gwas']['target_loci']

        # processed/default/trait/Whole_Blood/EUR/sqtl/
        self.tool_parent_dir = os.path.join(self.study_dir, 
                                            self.gwas_trait, 
                                            self.biological_context, 
                                            self.population,
                                            self.qtl_type)
        
        Path(self.tool_parent_dir).mkdir(exist_ok=True, parents=True)

        self.min_matching_number = self.global_config['min_matching_number']

        self.rank_dir = os.path.join(self.tool_parent_dir, 'rank') 
        self.output_processed_dir
        Path(self.rank_dir).mkdir(exist_ok=True, parents=True)

        self.report_path = os.path.join(self.study_dir) # trait/processed/default/

        self.gwas_col_dict = self.global_config['input']['gwas']['col_name_mapping']
        self.qtl_col_dict = self.global_config['input']['qtl']['col_name_mapping']

        if self.gwas_col_dict['variant_id'] == None:
            self.gwas_col_dict['variant_id'] = 'variant_id'

        if ('qtl_preprocessed_dir' in self.global_config['input']['qtl']) and \
            self.global_config['input']['qtl']['qtl_preprocessed_dir']:
            # 1. use custom eqtl preprocessed file
            custom_preprocess_qlt_path = self.global_config['input']['qtl']['qtl_preprocessed_dir']

            # e.g. preprocessed/eqtl/Whole_Blood/grouped
            self.qtl_grouped_dir = os.path.join(custom_preprocess_qlt_path, 
                                                self.qtl_type, 
                                                self.biological_context, 
                                                'grouped') 
            
            # e.g. preprocessed/eqtl/Whole_Blood/filtered_gene.tsv.gz
            self.qtl_output_report = os.path.join(custom_preprocess_qlt_path, 
                                                   self.qtl_type, 
                                                   self.biological_context, 
                                                   'filtered_gene.tsv.gz')
            self.qtl_LD_window = os.path.join(custom_preprocess_qlt_path, 
                                                   self.qtl_type, 
                                                   self.biological_context, 
                                                   'LD_window_per_gene.tsv.gz')
            self.qtl_preprocesed_dir = os.path.join(custom_preprocess_qlt_path,
                                                    self.qtl_type, 
                                                   self.biological_context)
        else:
            # 2. set eqtl output path
            self.qtl_grouped_dir = os.path.join(self.output_preprocessed_dir, 
                                                self.qtl_type, 
                                                self.biological_context, 
                                                'grouped')
            self.qtl_output_report = os.path.join(self.output_preprocessed_dir, 
                                                   self.qtl_type, 
                                                   self.biological_context,
                                                   'filtered_gene.tsv.gz')
            self.qtl_LD_window = os.path.join(self.output_preprocessed_dir, 
                                              self.qtl_type, 
                                              self.biological_context, 
                                              'LD_window_per_gene.tsv.gz')
            self.qtl_preprocesed_dir = os.path.join(self.output_preprocessed_dir,
                                                    self.qtl_type, 
                                                   self.biological_context)

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
        print(f"self.ref_vcf_dir: {self.ref_vcf_dir}")

        self.parallel = parallel

        # p-val threshold
        self.gwas_p_threshold = self.global_config['p-value_threshold']['gwas']
        if self.gwas_p_threshold <= 0:
            self.gwas_p_threshold = 5.0E-8
        self.qtl_p_threshold = self.global_config['p-value_threshold']['qtl']
        if self.qtl_p_threshold <= 0:
            self.qtl_p_threshold = 1.0E-5
        # input sep
        self.gwas_sep = self.global_config['input']['gwas'].get('sep', '\t')
        if self.gwas_sep == '\\t':
            self.gwas_sep = '\t'
        elif (len(self.gwas_sep)) > 1:
            logging.warning(f'GWAS file separator is too long (can only be one char), using tab instead')
            self.gwas_sep = '\t'
        self.qtl_sep = self.global_config['input']['qtl'].get('sep', '\t')
        if self.qtl_sep == '\\t':
            self.qtl_sep = '\t'
        elif (len(self.qtl_sep)) > 1:
            logging.warning(f'QTL file separator is too long (can only be one char), using tab instead')
            self.qtl_sep = '\t'

        # tools parameter config file
        self.tools_config_file = tools_config_file
