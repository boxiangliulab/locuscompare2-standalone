import os
import gzip
import sqlite3
from pathlib import Path

import pandas as pd

SNP_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

weight_ddl = '''CREATE TABLE "weights" (
"gene" TEXT,
  "rsid" TEXT,
  "varID" TEXT,
  "ref_allele" TEXT,
  "eff_allele" TEXT,
  "weight" REAL
)'''

weight_idx1 = 'CREATE INDEX weights_gene ON weights (gene)'
weight_idx2 = 'CREATE INDEX weights_rsid ON weights (rsid)'
weight_idx3 = 'CREATE INDEX weights_rsid_gene ON weights (rsid, gene)'

extra_ddl = '''CREATE TABLE "extra" (
"gene" TEXT,
  "genename" TEXT,
  "gene_type" TEXT,
  "n.snps.in.model" INTEGER,
  "pred.perf.R2" TEXT,
  "pred.perf.pval" TEXT,
  "pred.perf.qval" TEXT
)'''


def create_genotype_from_vcf(input_vcf, output_dir):
    """
    Reference from split_genotype_by_chr.py
    Output files are output_dir/geno_chr{num}.tsv
    chrom in varID must be prefixed with 'chr'.
    """
    print('Create genotype files from vcf start')
    if not os.path.exists(input_vcf) or os.path.getsize(input_vcf) <= 0:
        raise ValueError(f'{input_vcf} does not exist or is empty')
    for line in gzip.open(input_vcf, 'rb'):
        line = line.decode('UTF-8').strip('\n')
        if line.startswith('#CHROM'):
            vcf_header = line.split('\t')
            break
    else:
        raise ValueError(f'{input_vcf} does not have a header row start with #CHROM')
    vcf_df = pd.read_table(input_vcf, header=None, comment='#')
    vcf_df.columns = vcf_header
    vcf_df.replace({'0/0': '0', '0/1': '1', '1/0': '1', '1/1': '2', './.': pd.NA}, inplace=True)
    vcf_chr_notation = vcf_df['#CHROM'].astype(str).str.contains('chr').any()
    # COUNTED column is major allele, .traw file
    vcf_df['varID'] = \
        vcf_df['#CHROM'].astype(str) if vcf_chr_notation else 'chr' + vcf_df['#CHROM'].astype(str) \
                                                              + '_' + vcf_df['POS'].astype(str) \
                                                              + '_' + vcf_df['REF'].astype(str) \
                                                              + '_' + vcf_df['ALT'].astype(str) \
                                                              + '_b38'
    indel_bool_series = (vcf_df['ALT'].str.len() != 1) | (vcf_df['REF'].str.len() != 1)
    vcf_df.drop(labels=vcf_df[indel_bool_series].index, inplace=True)
    ambiguous_strand_series = (vcf_df['ALT'] == vcf_df['REF'].map(SNP_COMPLEMENT))
    vcf_df.drop(labels=vcf_df[ambiguous_strand_series].index, inplace=True)
    vcf_df.drop_duplicates(subset='varID', inplace=True)
    if vcf_df.shape[0] == 0:
        raise ValueError(f'No eligible data in {input_vcf}')
    vcf_df.drop(columns=['POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
    grouped = vcf_df.groupby('#CHROM')
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for name, geno_group in grouped:
        geno_group.drop(columns='#CHROM', inplace=True)
        geno_group = geno_group.reindex(columns=['varID'] + [col for col in geno_group.columns if col not in ['varID']],
                                        copy=False)
        chr_notation = f'{name}' if vcf_chr_notation else f'chr{name}'
        output_file = os.path.join(output_dir, f'geno_{chr_notation}.tsv')
        geno_group.to_csv(output_file, sep='\t', header=True, index=False, na_rep='NA')
    print('Create genotype files from vcf complete')


def main(src_dir, dest_dir, genotype_dir, gene_info):
    gtf = pd.read_table(gene_info)
    for src_file in os.listdir(src_dir):
        if src_file.startswith('.') or os.path.isdir(src_file):
            continue
        if not src_file.endswith('.enet.wts'):
            continue
        gene_id = src_file.strip('.enet.wts')
        wts_path = os.path.join(src_dir, src_file)
        wts = []
        var_ids = []
        for line in open(wts_path):
            if not line.strip():
                continue
            wts_tuple = line.split()
            wts_tuple2 = wts_tuple[0].split('_')
            var_ids.append(wts_tuple[0])
            wts.append({"gene": gene_id,
                        "rsid": wts_tuple[0],
                        "varID": wts_tuple[0],
                        "ref_allele": wts_tuple2[2],
                        "eff_allele": wts_tuple2[3],
                        "weight": wts_tuple[1]})
        if len(wts) == 0:
            continue
        genotype_file = os.path.join(genotype_dir, f'{gene_id}.vcf.gz')
        cur_gene_info = gtf[gtf['gene_id'] == gene_id].iloc[0]
        extra = [{'gene': gene_id, 'genename': cur_gene_info.loc['gene_name'],
                  'gene_type': cur_gene_info.loc['gene_type'], 'n': len(wts)}]
        db_dir = os.path.join(dest_dir, gene_id, 'db')
        cov_dir = os.path.join(dest_dir, gene_id, 'covariances')
        Path(db_dir).mkdir(parents=True, exist_ok=True)
        Path(cov_dir).mkdir(parents=True, exist_ok=True)
        db_file = os.path.join(db_dir, 'result.db')
        cov_file = os.path.join(cov_dir, 'covariances.txt')
        if os.path.exists(db_file):
            os.remove(db_file)
        if os.path.exists(cov_file):
            os.remove(cov_file)
        with sqlite3.connect(db_file) as con:
            cur = con.cursor()
            cur.execute(weight_ddl)
            cur.executemany('INSERT INTO weights VALUES(:gene, :rsid, :varID, :ref_allele, :eff_allele, :weight)',
                            wts)
            con.commit()
            cur.execute(weight_idx1)
            cur.execute(weight_idx2)
            cur.execute(weight_idx3)

            cur.execute(extra_ddl)

            cur.executemany(
                'INSERT INTO extra(gene, genename, gene_type, "n.snps.in.model") VALUES (:gene, :genename, :gene_type, :n)',
                extra)
            con.commit()
        con.close()
        create_genotype_from_vcf(genotype_file, os.path.join(dest_dir, gene_id))
        genotype_dosage_file = os.path.join(dest_dir, gene_id, f'geno_chr{cur_gene_info.loc["chrom"]}.tsv')
        df2 = pd.read_table(genotype_dosage_file)
        df2 = df2[df2['varID'].isin(var_ids)]
        df2.set_index('varID', inplace=True)
        cov = df2.T.cov()
        with open(cov_file, 'w') as file:
            # Write the header
            file.write("GENE RSID1 RSID2 VALUE\n")

            # Write the upper triangle of the covariance matrix to the file
            columns = cov.columns
            for i in range(len(columns)):
                for j in range(i, len(columns)):
                    id1 = columns[i]
                    id2 = columns[j]
                    c = cov.iloc[i, j]
                    file.write(f"{gene_id} {id1} {id2} {c}\n")
        os.system(f'gzip {cov_file}')


if __name__ == '__main__':
    main('/Volumes/HD/biodata/sem/sim_scheme/sim_rst/predixcan_wts',
         '/Volumes/HD/biodata/sem/pred_model',
         '/Volumes/HD/biodata/sem/eqtl',
         '/Volumes/HD/biodata/colocalization-tools/raw/genecode/gtf_all_gene.tsv')
