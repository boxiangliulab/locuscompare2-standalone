import os
from pathlib import Path

import pandas as pd
import scipy.stats as ss
import sys
import ranking.intact as intact
import ranking.rra as rra


def run_ranking(rpt_obj=None, output_file_path=None, prior_fun=None, sample_size=None):
    
    
    out_dir = os.path.dirname(output_file_path)
    if not os.path.exists(out_dir):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
    geo_output_file = os.path.join(out_dir, 'geo.tsv')
    intact_output_file = os.path.join(out_dir, 'intact.tsv')
    geo_result = rra.run_ranking(output_file_path=geo_output_file, rpt=rpt_obj, sample_size=sample_size, method='GEO')
    intact_result = intact.run_ranking(rpt_obj, intact_output_file, prior_fun)

    geo_result_exist = geo_result is not None and os.path.exists(geo_result) and os.path.getsize(geo_result) > 0
    intact_result_exist = intact_result is not None and os.path.exists(intact_result) and os.path.getsize(
        intact_result) > 0

    if geo_result_exist and intact_result_exist:
        # output cols: ['gene_id', 'geo_ranking', 'avg_ranking', 'intact_probability']
        rgeo_df = pd.read_table(geo_result, usecols=['gene_id', 'geo_ranking'])
        intact_df = pd.read_table(intact_result,
                                  usecols=lambda col: col in ['gene_id', 'avg_ranking', 'intact_probability'])
        result_df = pd.merge(left=rgeo_df, right=intact_df,
                             left_on='gene_id', right_on='gene_id',
                             how='outer')
        result_df['geo_ranking'] = ss.rankdata(result_df['geo_ranking'])
        result_df['geo_ranking'] = result_df['geo_ranking'].apply(lambda x: int(x))
        result_df['avg_ranking'] = ss.rankdata(result_df['avg_ranking'])
        result_df['avg_ranking'] = result_df['avg_ranking'].apply(lambda x: int(x))
    elif geo_result_exist:
        # intact_probability column is not present if user run TWAS only or colocalization methods
        result_df = pd.read_table(geo_result, usecols=['gene_id', 'geo_ranking'])
        result_df['geo_ranking'] = ss.rankdata(result_df['geo_ranking'])
        result_df['geo_ranking'] = result_df['geo_ranking'].apply(lambda x: int(x))
        os.remove(geo_result)
    elif intact_result_exist:
        result_df = pd.read_table(intact_result,
                                  usecols=lambda col: col in ['gene_id', 'avg_ranking', 'intact_probability'])
        result_df['avg_ranking'] = ss.rankdata(result_df['avg_ranking'])
        result_df['avg_ranking'] = result_df['avg_ranking'].apply(lambda x: int(x))
        os.remove(intact_result)
    else:
        result_df = None

    if result_df is not None:
        result_df.to_csv(output_file_path, sep='\t', header=True, index=False, na_rep='NA')

    return output_file_path
