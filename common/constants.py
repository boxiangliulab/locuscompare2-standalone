from pathlib import Path
import os

column_spliter = '\t'
output_spliter = '\t'
gene_id_prefix = 'ENSG'
manhattan_min_pval = 0.5
manhattan_max_plot = 50000
default_config = os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), 'config.yml')
default_study = 'default'
default_report_gene_ranking = 50
# tools_config = os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), 'tools_config.yml')
SNP_ALLELE = ['A', 'T', 'C', 'G']
