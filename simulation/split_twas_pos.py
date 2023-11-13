import pandas as pd

if __name__ == '__main__':
    df = pd.read_table('/Volumes/HD/biodata/sem/twas_model/test_tissue.pos')
    for i in range(0, df.shape[0]):
        gene_id = df.at[i, 'ID']
        df.loc[[i], ].to_csv(f'{gene_id}.pos', index=False, header=True, sep='\t')
