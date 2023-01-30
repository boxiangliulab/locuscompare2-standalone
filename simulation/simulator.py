import argparse
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd

import pheno_to_bed


# generate output_prefix.samples, output_prefix.hap, output_prefix.legend, output_prefix_maf.tsv
# in current dir
def prepare_hapgen2_input(subset_vcf, output_prefix, maf_thresh):
    print(f'converting subset_vcf to haplotype files')
    # generate output_prefix.samples, output_prefix.hap.gz, output_prefix.legend.gz
    # bcftools convert 和 bcftools view对于-i 'TYPE="snp"'的执行结果不同, view结果包含MULTI_ALLELIC
    os.system(f'bcftools convert -i \'TYPE="snp"\' -h {output_prefix} --vcf-ids {subset_vcf}')
    # generate output_prefix.hap, output_prefix.legend
    os.system(f'gunzip -f {output_prefix}.hap.gz {output_prefix}.legend.gz')
    # generate output_prefix_pre.frq, this will be deleted
    print(f'calculating MAF from vcf')
    os.system(f'plink --vcf {subset_vcf} --freq --out {output_prefix}_pre')
    maf_df = pd.read_table(f'{output_prefix}_pre.frq', sep=r'\s+')
    # include all variant that has maf>maf_thresh
    maf_df = maf_df[maf_df['MAF'] > maf_thresh]
    # vcf ID column must NOT be missing
    vcf_df = pd.read_table(subset_vcf, header=None, comment='#', usecols=[1, 2, 3, 4])
    vcf_df.columns = ['POS', 'SNP', 'REF', 'ALT']
    maf_df = pd.merge(maf_df, vcf_df,
                      left_on='SNP',
                      right_on='SNP',
                      how='left')
    # choose which alt===A1 1 else 0
    maf_df['ALLELE_STATE'] = (maf_df['ALT'] == maf_df['A1']).map({True: 1, False: 0})
    maf_df.drop(['REF', 'ALT'], axis=1, inplace=True)
    maf_df.drop(labels=maf_df[maf_df['SNP'] == '.'].index, inplace=True)
    if maf_df.shape[0] < 2:
        print(f'warning: maf_df should has at least 2 rows')
        return None
    maf_df.to_csv(f'{output_prefix}_maf.tsv', sep='\t', header=True, index=False)
    return output_prefix, f'{output_prefix}_maf.tsv'


# generate control data in plink binary format: output_prefix_ctrl.bed, output_prefix_ctrl.bim, output_prefix_ctrl.fam
# case data in plink binary format: output_prefix_case.bed, output_prefix_case.bim, output_prefix_case.fam
# merged case/control data in binary format: output_prefix.bed, output_prefix.bim, output_prefix.fam
# allele encoding will be converted from hapgen2 format to plink format: 0/1 to 1/2
def generate_cc(genetic_map_path, input_prefix, chrom, disease_loci_pos_and_allele_state_list, ctrl_cnt, case_cnt,
                output_prefix, merge_cc=True, disease_odds_ratio=1.1, effect_population_size=11418):
    disease_param = ''
    print(f'generating case control data, disease loci are: {disease_loci_pos_and_allele_state_list}')
    for disease_pos_and_allele_state in disease_loci_pos_and_allele_state_list:
        disease_param += f' {disease_pos_and_allele_state} {disease_odds_ratio} {disease_odds_ratio * disease_odds_ratio}'
    os.system(f'hapgen2 -m {genetic_map_path} '
              f'-l {input_prefix}.legend '
              f'-h {input_prefix}.hap '
              f'-o {output_prefix} '
              f'-dl {disease_param} '
              f'-n {ctrl_cnt} {case_cnt} '
              # f'-output_snp_summary '
              f'-Ne {effect_population_size}')
    print(f'converting case control data to plink binary format')
    os.system(f'plink --data {output_prefix}.controls '
              f'-oxford-single-chr {chrom} --allow-no-sex --make-bed --out {output_prefix}_ctrl')
    os.system(f'plink --data {output_prefix}.cases '
              f'-oxford-single-chr {chrom} --allow-no-sex --make-bed --out {output_prefix}_case')
    if merge_cc:
        os.system(f'plink --bfile {output_prefix}_ctrl --bmerge {output_prefix}_case --out {output_prefix}')
        # delete intermediate case/control binary files
        os.remove(f'{output_prefix}_ctrl.bed')
        os.remove(f'{output_prefix}_ctrl.bim')
        os.remove(f'{output_prefix}_ctrl.fam')
        os.remove(f'{output_prefix}_case.bed')
        os.remove(f'{output_prefix}_case.bim')
        os.remove(f'{output_prefix}_case.fam')
    # delete intermediate hapgen2 output case/control files
    os.remove(f'{output_prefix}.controls.sample')
    os.remove(f'{output_prefix}.controls.haps')
    os.remove(f'{output_prefix}.controls.gen')
    os.remove(f'{output_prefix}.cases.sample')
    os.remove(f'{output_prefix}.cases.haps')
    os.remove(f'{output_prefix}.cases.gen')
    return output_prefix if merge_cc else f'{output_prefix}_case', f'{output_prefix}_ctrl'


# generate output_prefix.frq, output_prefix.assoc, output_prefix.assoc.logistic, output_prefix.log
def perform_gwas_analysis(input_prefix, output_id, maf_thresh=0.05):
    print(f'performing gwas analysis')
    # QC steps:
    # Investigate missingness per individual and per SNP
    # Delete SNPs with missingness >0.02
    # Delete individuals with missingness >0.02
    # SKIP: Check sex discrepancy
    # MAF check
    os.system(f'plink --allow-no-sex '
              f'--bfile {input_prefix} '
              f'--geno 0.02 '
              f'--mind 0.02 '
              f'--maf {maf_thresh} --make-bed '
              f'--out {output_id}_qc')
    # SKIP: HWE check
    # SKIP: Heterozygosity rate check
    # SKIP: Cryptic relatedness check
    # SKIP: Population stracfication check?
    # GWAS analysis
    os.system(f'plink --allow-no-sex '
              f'--bfile {output_id}_qc '
              f'--freq --logistic beta --assoc --ci 0.95 '
              f'--out {output_id}')
    # regression was tested against A1, i.e. minor allele
    # since effect allele is the allele whose effects in relation to disease are being studied
    # so here A1=minor allele=effect allele, EAF=MAF
    gwas_maf_df = pd.read_table(f'{output_id}.frq', sep=r'\s+', usecols=['SNP', 'A2', 'MAF'])
    gwas_df = pd.read_table(f'{output_id}.assoc.logistic', header=0, sep=r'\s+')
    gwas_df = pd.merge(gwas_df, gwas_maf_df,
                       left_on='SNP',
                       right_on='SNP',
                       how='left')
    gwas_df.to_csv(f'gwas_{output_id}.tsv', sep='\t', header=True, index=False)
    return f'gwas_{output_id}.tsv'


# generate eqtl individual data: output_bed_prefix.bed.gz, output_vcf_prefix.vcf.gz
def generate_eqtl(input_prefix, causal_snplist_path, output_bed_prefix,
                  output_vcf_prefix, chrom, start_pos, end_pos,
                  gene_id, causal_mode, strand=None, heritability=0.2, simu_rep=1):
    print(f'generating eqtl expression data')
    # generate phenotype data
    pheno_output_prefix = os.path.join(os.path.dirname(output_bed_prefix),
                                       f'{output_bed_prefix.strip("_exp")}_gcta_pheno')
    os.system(f'gcta --bfile {input_prefix}  '
              f'--simu-qt --simu-hsq {0 if causal_mode == 0 else heritability} '
              f'--simu-causal-loci {causal_snplist_path} --simu-rep {simu_rep} '
              f'--out {pheno_output_prefix}')
    # convert pheno data to bed format and then index it
    pheno_to_bed.convert_gcta_pheno_to_bed(f'{pheno_output_prefix}.phen', f'{output_bed_prefix}.bed',
                                           chrom, start_pos, end_pos, gene_id, strand)
    os.system(f'bgzip -f {output_bed_prefix}.bed && tabix -f -p bed {output_bed_prefix}.bed.gz')
    # convert input genotype data to vcf format, use original uncompressed vcf to correct allele order and then index it
    os.system(f'plink --allow-no-sex --bfile {input_prefix} '
              f'--real-ref-alleles '
              f'--a2-allele {gene_id}.vcf 4 3 \'#\' '
              f'--recode vcf-iid tab bgz '
              f'--out {output_vcf_prefix}')
    os.system(f'tabix -f -p vcf {output_vcf_prefix}.vcf.gz')
    return f'{output_vcf_prefix}.vcf.gz', f'{output_bed_prefix}.bed.gz'


def perform_eqtl_analysis(phenotype_bed_path, genotype_vcf_path, output_path, p_val_thresh=0.1):
    print(f'performing eqtl analysis')
    # use first 3 PC as cov
    print(f'performing pca analysis')
    # gcta --simu-rep几次就会有几行pca结果
    chunk_count = get_bed_chrom_count(phenotype_bed_path)
    cov_param = ''
    if chunk_count > 1:
        os.system(f'QTLtools pca --bed {phenotype_bed_path} --out phenotype_raw_pca --center --scale && '
                  f'head -4 phenotype_raw_pca.pca > phenotype_pca.pca')
        cov_param = '--cov phenotype_pca.pca '
    for chunk in range(0, chunk_count + 1):
        os.system(f'QTLtools cis '
                  f'--vcf {genotype_vcf_path} '
                  f'--bed {phenotype_bed_path} '
                  f'{cov_param}'
                  f'--nominal {p_val_thresh} '
                  f'--normal --std-err '
                  f'--chunk {chunk} {chunk_count} '
                  f'--out qtl_result_{chunk}')
    qtl_result_list = []
    for f in os.listdir(os.getcwd()):
        if f == 'qtl_result_0':
            continue
        if f.startswith('qtl_result_'):
            qtl_result_list.append(pd.read_table(f, header=None, sep=r'\s+'))
    eqtl_df = pd.concat(qtl_result_list)
    # different mode has different result column
    eqtl_df.columns = ['phe_id', 'phe_chrom', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var',
                       'var_id', 'var_chrom', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'beta', 'se', 'best_hit']
    # calculate MAF for eqtl data
    os.system(f'plink --vcf {genotype_vcf_path} --freq --double-id --out  eqtl_maf')
    eqtl_maf_df = pd.read_table(f'eqtl_maf.frq', sep=r'\s+', usecols=['SNP', 'MAF'])
    eqtl_maf_df.columns = ['var_id', 'maf']
    eqtl_df = pd.merge(eqtl_df, eqtl_maf_df,
                       left_on='var_id',
                       right_on='var_id',
                       how='left')
    # merge REF, ALT from genotype vcf
    genotype_df = pd.read_table(genotype_vcf_path, header=None, comment='#', usecols=[2, 3, 4])
    genotype_df.columns = ['var_id', 'ref', 'alt']
    eqtl_df = pd.merge(eqtl_df, genotype_df,
                       left_on='var_id',
                       right_on='var_id',
                       how='left')
    eqtl_df.to_csv(output_path, sep='\t', header=True, index=False)
    return output_path


def simulate_for_chrom(input_vcf, chrom_num_in_vcf, gwas_causal_list, eqtl_causal_list, gene_list_path,
                       generated_file_list, genetic_map, ctrl_count, case_count, maf_thresh,
                       eqtl_genetic_models, heritability=0.2, disease_risk=1.1):
    start_time = datetime.now()
    print(f'simulate for chrom {chrom_num_in_vcf} start at {start_time}')
    gene_list_df = pd.read_table(gene_list_path)
    gene_list_df = gene_list_df.loc[gene_list_df['chrom'] == chrom_num_in_vcf]
    print(f'gene count for chromosome {chrom_num_in_vcf} is: {gene_list_df.shape[0]}')
    # sampling ranges could overlap with each other
    # gene average len: 31kb
    gene_len_thresh = 200000
    for _, row in gene_list_df.iterrows():
        gene_len = row.loc['end'] - row.loc['start']
        distance_to_include = round((gene_len_thresh - gene_len) / 2) if gene_len < gene_len_thresh else 0
        start = max(row.loc['start'] - distance_to_include, 1)
        end = row.loc['end'] + distance_to_include
        chrom = row.loc['chrom']
        gene_id = row.loc['gene_id']
        center_snp = row.loc['center_snp']
        # keep only snps from vcf
        print(f'extracting and filtering subset vcf for gene {gene_id}')
        chrom_code_in_vcf = f'chr{chrom}' if is_vcf_chrom_code_contains_chr(input_vcf) else chrom
        pre_filtered_vcf = f'{gene_id}_pre_filtered.vcf.gz'
        with open(pre_filtered_vcf, mode='w') as filtered_vcf:
            proc1 = subprocess.Popen(['tabix', '-h', input_vcf, f'{chrom_code_in_vcf}:{start}-{end}'],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc2 = subprocess.Popen(['bcftools', 'view', '-m2', '-M2', '-v', 'snps', '-Ou'],
                                     stdin=proc1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc1.stdout.close()  # Allow pro1 to receive a SIGPIPE if proc2 exits.
            proc3 = subprocess.Popen(['bcftools', 'norm', '-d', 'any', '-Oz'],
                                     stdin=proc2.stdout, stdout=filtered_vcf, stderr=subprocess.PIPE)
            proc2.stdout.close()  # Allow proc2 to receive a SIGPIPE if proc3 exits.
            # bcftools output summary numbers to stderr,
            # so can not use stderr of proc3.communicate to check if subproces is success
            # https://github.com/samtools/bcftools/issues/850
            proc1.wait()
            proc2.wait()
            proc3.communicate()
        os.system(f'tabix -f -p vcf {pre_filtered_vcf}')
        if get_vcf_records_count(pre_filtered_vcf) == 0:
            print(f'warning: no variant in {pre_filtered_vcf}')
            continue
        # set varid for variant that still missing rsid
        # keep *annotated* *uncompressed* vcf for later allele order ref in generating qtltools genotype vcf
        os.system(f'bcftools annotate -r {chrom_code_in_vcf} -I +"%CHROM:%POS\\_%REF\\_%FIRST_ALT" '
                  f'{pre_filtered_vcf} -Ov -o {gene_id}.vcf && bgzip -f -k {gene_id}.vcf')
        os.remove(pre_filtered_vcf)
        os.remove(f'{pre_filtered_vcf}.tbi')
        print(f'set var_id for missing SNP completed')
        # os.remove(f'{gene_id}_anno1.vcf')
        # generate hapgen2 required files
        haps_leg_maf_tuple = prepare_hapgen2_input(f'{gene_id}.vcf.gz', gene_id, maf_thresh)
        if haps_leg_maf_tuple is None:
            continue
        disease_loci_pos_and_allele_state_list = []
        maf_df = pd.read_table(f'{gene_id}_maf.tsv')
        # pick the center variant as disease loci
        # disease_loci is a series
        disease_loci = maf_df[maf_df['SNP'] == center_snp].iloc[0]
        with open(gwas_causal_list, mode='a') as gwas_causal:
            gwas_causal.write(f'{chrom}\t{gene_id}\t{disease_loci.loc["SNP"]}\t{disease_loci.loc["POS"]}\n')
        disease_loci_pos_and_allele_state_list.append(f'{disease_loci.loc["POS"]} {disease_loci.loc["ALLELE_STATE"]}')
        # use hapgen2 to generate case control data
        gwas_simulation_prefix = f'{gene_id}_gwas'
        generate_cc(genetic_map,
                    gene_id,
                    chrom,
                    disease_loci_pos_and_allele_state_list,
                    ctrl_count, case_count,
                    gwas_simulation_prefix,
                    disease_odds_ratio=disease_risk)
        # by default, plink will filter out records with ld r2<0.2, --ld-window-r2 0 prevent this behaviour
        # monomorphic variant LD r2 will not be present in report in table(inter-chr) format,
        # in square format output they will be nan
        os.system(f'plink --vcf {gene_id}.vcf.gz '
                  f'--r2 inter-chr with-freqs --ld-window-r2 0 '
                  f'--ld-snp {disease_loci.loc["SNP"]} '
                  f'--out {gene_id}')
        # the above cmd outputs table with cols: CHR_A,BP_A,SNP_A,MAF_A,CHR_B,BP_B,SNP_B,MAF_B,R2
        # in which CHR_A = disease loci snp
        if not Path(f'{gene_id}.ld').exists():
            print(f'warning: no ld generated for {gene_id}.vcf.gz, are the variants in vcf monomorphic?')
            continue
        ld_df = pd.read_table(f'{gene_id}.ld', sep=r'\s+', header=0)
        sample_gwas_summary_stats = perform_gwas_analysis(gwas_simulation_prefix, f'{gene_id}')
        gwas_df = pd.read_table(sample_gwas_summary_stats)
        gwas_df.sort_values(['P'], inplace=True)
        gwas_max_assoc_row = gwas_df.iloc[0]
        if center_snp == gwas_max_assoc_row.loc["P"]:
            gwas_causal_max_assoc_r2 = 1
        else:
            gwas_causal_max_assoc_r2 = ld_df[ld_df['SNP_B'] == gwas_max_assoc_row["SNP"]].iloc[0].loc['R2']
        # drop the row which represent r2 with self(r2=1)
        ld_df.drop(index=ld_df.loc[ld_df['SNP_B'] == ld_df['SNP_A']].index, inplace=True)
        # keep variants with maf > maf_thresh
        ld_df.drop(index=ld_df.loc[ld_df['MAF_B'] <= maf_thresh].index, inplace=True)
        # keep variants within 50kb of disease loci
        eqtl_gwas_causal_distance = 50000
        ld_df.drop(index=ld_df.loc[(ld_df['BP_A'] - ld_df['BP_B']).abs() > eqtl_gwas_causal_distance].index,
                   inplace=True)
        generated_file_row = f'{chrom}\t{gene_id}\t{start}\t{end}\t{center_snp}'
        generated_file_row += f'\t{gwas_simulation_prefix}\t{sample_gwas_summary_stats}\t{gwas_max_assoc_row["SNP"]}'
        generated_file_row += f'\t{gwas_max_assoc_row.loc["P"]}\t{gwas_causal_max_assoc_r2}'
        for eqtl_causal_type in eqtl_genetic_models:
            # generate 250 control data for eQTL
            eqtl_simulation_prefix = f'{gene_id}_eqtl_{eqtl_causal_type}'
            generate_cc(genetic_map,
                        gene_id,
                        chrom,
                        disease_loci_pos_and_allele_state_list,
                        250, 1,
                        eqtl_simulation_prefix,
                        False,
                        disease_odds_ratio=1)
            # delete case binary data generated in hapgen2 eqtl individual genotype simulation
            os.remove(f'{eqtl_simulation_prefix}_case.bed')
            os.remove(f'{eqtl_simulation_prefix}_case.bim')
            os.remove(f'{eqtl_simulation_prefix}_case.fam')
            # eqtl_loci_causal_list_file是传给gcta的当前模拟区段的causal, eqtl_causal_list是集合了所有causal list
            eqtl_loci_causal_list_file = f'{gene_id}_eqtl_causal_{eqtl_causal_type}_list.txt'
            with open(eqtl_loci_causal_list_file, mode='w') as loci_causal_ptr, \
                    open(eqtl_causal_list, mode='a') as eqtl_causal:
                if eqtl_causal_type in [0, 1]:
                    # no causal variant or same causal variant as gwas
                    eqtl_causal.write(
                        f'{chrom}\t{gene_id}\t{disease_loci.loc["SNP"]}\t{disease_loci.loc["POS"]}\t{eqtl_causal_type}\n')
                    loci_causal_ptr.write(f'{disease_loci.loc["SNP"]}\n')
                else:
                    # different causal variant as gwas
                    # keep those variants that: r2_min < r2 <= r2_max
                    if eqtl_causal_type == 2:
                        r2_min = 0
                        r2_max = 0.4
                    elif eqtl_causal_type == 3:
                        r2_min = 0.4
                        r2_max = 0.7
                    else:
                        r2_min = 0.7
                        r2_max = 0.9
                    sub_ld_df = ld_df[(ld_df['R2'] > r2_min) & (ld_df['R2'] <= r2_max)]
                    if sub_ld_df.shape[0] == 0:
                        print(f'warning: no suitable eQTL causal variant with MAF>{maf_thresh} and '
                              f'range in {(r2_min, r2_max)} of LD r2 within 50kb of disease loci {disease_loci.loc["SNP"]}')
                        generated_file_row += f'\tNA\tNA\tNA\tNA\tNA'
                        continue
                    # for H2, a single cohort with up to five distinct causal variant
                    if sub_ld_df.shape[0] < 5:
                        eqtl_causal_loci_df = sub_ld_df
                    else:
                        eqtl_causal_loci_df = sub_ld_df.sample(5)
                    for _, eqtl_loci_row in eqtl_causal_loci_df.iterrows():
                        eqtl_causal.write(
                            f'{chrom}\t{gene_id}\t{eqtl_loci_row.loc["SNP_B"]}\t{eqtl_loci_row.loc["BP_B"]}\t{eqtl_causal_type}\n')
                        loci_causal_ptr.write(f'{eqtl_loci_row.loc["SNP_B"]}\n')
            geno_vcf, pheno_bed = generate_eqtl(f'{eqtl_simulation_prefix}_ctrl', eqtl_loci_causal_list_file,
                                                f'{gene_id}_exp_{eqtl_causal_type}',
                                                f'{gene_id}_geno_{eqtl_causal_type}',
                                                chrom, row.loc['start'], row.loc['end'], gene_id, eqtl_causal_type,
                                                row.loc['strand'], heritability=heritability)
            sample_eqtl_summary_stats = perform_eqtl_analysis(pheno_bed, geno_vcf,
                                                              f'eqtl_{eqtl_causal_type}_{gene_id}.tsv')
            eqtl_df = pd.read_table(sample_eqtl_summary_stats)
            eqtl_df.sort_values(['nom_pval'], inplace=True)
            eqtl_max_assoc_row = eqtl_df.iloc[0]
            generated_file_row += f'\t{gene_id}_geno_{eqtl_causal_type}.vcf.gz\t{gene_id}_exp_{eqtl_causal_type}.bed.gz'
            generated_file_row += f'\t{sample_eqtl_summary_stats}\t{eqtl_max_assoc_row.loc["var_id"]}'
            generated_file_row += f'\t{eqtl_max_assoc_row.loc["nom_pval"]}'
            # delete control binary data generated in hapgen2 eqtl individual genotype simulation
            os.remove(f'{eqtl_simulation_prefix}_ctrl.bed')
            os.remove(f'{eqtl_simulation_prefix}_ctrl.bim')
            os.remove(f'{eqtl_simulation_prefix}_ctrl.fam')
        # delete hapgen2 input files
        os.remove(f'{gene_id}.legend')
        os.remove(f'{gene_id}.hap')
        os.remove(f'{gene_id}.samples')
        generated_file_row += '\n'
        with open(generated_file_list, mode='a') as generated_file:
            generated_file.write(generated_file_row)


def get_vcf_records_count(vcf_file):
    if not (Path(f'{vcf_file}.tbi').exists() or Path(f'{vcf_file}.csi').exists()):
        os.system(f'tabix -f -p vcf {vcf_file}')
    process = subprocess.Popen(['bcftools', 'index', '-n', vcf_file],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        return 0
    return int(stdout)


def is_vcf_chrom_code_contains_chr(vcf_file):
    if not (Path(f'{vcf_file}.tbi').exists() or Path(f'{vcf_file}.csi').exists()):
        os.system(f'tabix -f -p vcf {vcf_file}')
    process = subprocess.Popen(['tabix', '-l', vcf_file],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        return False
    return 'chr' in stdout.decode("UTF-8")


def get_bed_chrom_count(bed_file):
    if not (Path(f'{bed_file}.tbi').exists() or Path(f'{bed_file}.csi').exists()):
        os.system(f'tabix -f -p bed {bed_file}')
    process = subprocess.Popen(['tabix', '-l', bed_file],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        return False
    count = len(stdout.decode("UTF-8").splitlines())
    return count if count > 0 else 1


def simulate(output_suffix, src_vcf_dir, gene_list_file, genetic_map, ctrl_count, case_count, maf_thresh,
             heritability=0.2, disease_risk=1.1, eqtl_genetic_model=2, output_dir=None):
    start_time = datetime.now()
    print(f'Simulate start at: {start_time}')
    # files used to track causal list
    gwas_causal_list = f'gwas_causal_list_{output_suffix}.tsv'
    eqtl_causal_list = f'eqtl_causal_list_{output_suffix}.tsv'
    # files used to track gwas & eqtl geno/pheno files in generated order
    generated_file_list = f'generated_{output_suffix}.tsv'
    # write causal list header
    with open(gwas_causal_list, mode='w') as gwas_causal:
        gwas_causal.write(f'chrom\tgene\tsnp\tposition\n')
    with open(eqtl_causal_list, mode='w') as eqtl_causal:
        eqtl_causal.write(f'chrom\tgene\tsnp\tposition\tcausal_type\n')
    # write files list header
    eqtl_genetic_models = [1, eqtl_genetic_model]
    with open(generated_file_list, mode='w') as generated_file:
        header = f'chrom\tgene\tstart\tend\tgwas_causal_snp'
        header += f'\tgwas_geno_prefix\tgwas_sum_stats\tgwas_max_assoc_snp\tgwas_max_assoc_p\tgwas_r2'
        for i in eqtl_genetic_models:
            header += f'\teqtl_geno_{i}\teqtl_pheno_{i}\teqtl_sum_stats_{i}\teqtl_max_assoc_snp_{i}\teqtl_max_assoc_p_{i}'
        header += '\n'
        generated_file.write(header)

    for file in sorted(os.listdir(src_vcf_dir),
                       key=lambda f: int(re.search(r'\d+', f).group(0) if re.search(r'\d+', f) else 0)):
        if re.search(r'^chr\d+\.vcf\.gz$', file):
            chrom = int(re.search(r'\d+', file).group(0))
            simulate_for_chrom(os.path.join(src_vcf_dir, file), chrom, gwas_causal_list, eqtl_causal_list,
                               gene_list_file, generated_file_list, genetic_map, ctrl_count, case_count, maf_thresh,
                               eqtl_genetic_models, heritability=heritability, disease_risk=disease_risk)
    if output_dir is None:
        output_dir = f'output_{output_suffix}'
    Path(output_dir).mkdir(exist_ok=True, parents=True)
    shutil.copy(gene_list_file, output_dir)
    shutil.move(gwas_causal_list, output_dir)
    shutil.move(eqtl_causal_list, output_dir)
    shutil.move(generated_file_list, output_dir)
    merge_result = merge_sum_stat_files(os.path.join(output_dir, generated_file_list), output_suffix,
                                        eqtl_genetic_models=eqtl_genetic_models, output_dir=output_dir,
                                        heritability=heritability, disease_risk=disease_risk)
    merge_raw_files(os.path.join(output_dir, generated_file_list), output_suffix,
                    eqtl_genetic_models=eqtl_genetic_models, output_dir=output_dir)
    print(f'Simulate completed, duration {datetime.now() - start_time}')
    return merge_result


def merge_raw_files(generated_file_list, output_suffix, eqtl_genetic_models, output_dir):
    start_time = datetime.now()
    print(f'Merge raw files start at: {start_time}')
    file_list_df = pd.read_table(generated_file_list)
    file_list_df.sort_values(['chrom', 'start'], inplace=True)
    # concat eqtl genotype files
    for i in eqtl_genetic_models:
        os.system(
            f'bcftools concat -a -D {" ".join(file_list_df[f"eqtl_geno_{i}"].dropna())} -Oz -o eqtl_geno_{i}_{output_suffix}.vcf.gz '
            f'&& tabix -f -p vcf eqtl_geno_{i}_{output_suffix}.vcf.gz')
        # concat and write eqtl phenotype files
        eqtl_exp_list = []
        for indx, file_row in file_list_df.iterrows():
            if pd.isna(file_row.loc[f'eqtl_pheno_{i}']):
                continue
            eqtl_exp_list.append(pd.read_table(file_row.loc[f'eqtl_pheno_{i}'], sep=r'\s+'))
        if len(eqtl_exp_list) == 0:
            continue
        eqtl_exp_df = pd.concat(eqtl_exp_list)
        eqtl_exp_df.sort_values(['#chr', 'start'], inplace=True)
        eqtl_exp_df.to_csv(f'eqtl_exp_{i}_{output_suffix}.bed', sep='\t', header=True, index=False)
        os.system(f'bgzip -f eqtl_exp_{i}_{output_suffix}.bed && tabix -f -p bed eqtl_exp_{i}_{output_suffix}.bed.gz')
    # concat and write gwas case control files
    gwas_cc_binary = f'gwas_cc_{output_suffix}'
    gwas_list_df = file_list_df['gwas_geno_prefix'].copy()
    gwas_list_df.dropna(inplace=True)
    gwas_list_df.reset_index(inplace=True, drop=True)
    gwas_first_file = gwas_list_df.iloc[0]
    gwas_file_list = gwas_list_df.iloc[1:]
    gwas_file_list.to_csv(f'gwas_binary_to_merge.txt', sep='\t', header=False, index=False)
    os.system(f'plink --bfile {gwas_first_file} --merge-list gwas_binary_to_merge.txt --out {gwas_cc_binary}')
    print(f'Simulate completed, duration {datetime.now() - start_time}, output files are:')
    if output_dir is None:
        output_dir = f'output_{output_suffix}'
    Path(output_dir).mkdir(exist_ok=True, parents=True)
    for i in eqtl_genetic_models:
        shutil.move(f'eqtl_geno_{i}_{output_suffix}.vcf.gz', output_dir)
        shutil.move(f'eqtl_geno_{i}_{output_suffix}.vcf.gz.tbi', output_dir)
        shutil.move(f'eqtl_exp_{i}_{output_suffix}.bed.gz', output_dir)
        shutil.move(f'eqtl_exp_{i}_{output_suffix}.bed.gz.tbi', output_dir)
    shutil.move(f'{gwas_cc_binary}.bed', output_dir)
    shutil.move(f'{gwas_cc_binary}.bim', output_dir)
    shutil.move(f'{gwas_cc_binary}.fam', output_dir)
    print(f'Merge raw completed, duration {datetime.now() - start_time}')
    return (os.path.join(output_dir, gwas_cc_binary),
            [(os.path.join(output_dir, f'eqtl_geno_{i}_{output_suffix}.vcf.gz'),
              os.path.join(output_dir, f'eqtl_exp_{i}_{output_suffix}.bed.gz')) for i in eqtl_genetic_models])


def merge_sum_stat_files(generated_file_list, output_suffix, eqtl_genetic_models,
                         output_dir, heritability=0.2, disease_risk=1.1):
    start_time = datetime.now()
    print(f'Merge sumstat files start at: {start_time}')
    file_list_df = pd.read_table(generated_file_list)
    file_list_df.sort_values(['chrom', 'start'], inplace=True)
    # concat eqtl genotype files
    for i in eqtl_genetic_models:
        eqtl_sum_stats_list = []
        for indx, file_row in file_list_df.iterrows():
            if pd.isna(file_row.loc[f'eqtl_sum_stats_{i}']):
                continue
            if file_row.loc[f'eqtl_max_assoc_p_{i}'] > 0.01:
                # reject cohorts showing weak maximum association signals (P > 0.01 for eQTL cohorts)
                continue
            eqtl_sum_stats_list.append(pd.read_table(file_row.loc[f'eqtl_sum_stats_{i}'], sep=r'\s+'))
        if len(eqtl_sum_stats_list) == 0:
            print(f'warning: no suitable eQTL summary statistics file generated for causal type {i}')
            continue
        eqtl_sum_stats_df = pd.concat(eqtl_sum_stats_list)
        eqtl_sum_stats_df.sort_values('nom_pval', inplace=True)
        eqtl_sum_stats_df.to_csv(f'eqtl_simu_{i}_hsq{heritability}_{output_suffix}.tsv', sep='\t', header=True,
                                 index=False)
    gwas_sum_stats_list = []
    for indx, file_row in file_list_df.iterrows():
        if pd.isna(file_row.loc[f'gwas_sum_stats']):
            continue
        if file_row.loc[f'gwas_max_assoc_p'] > 1.0E-5:
            # reject cohorts showing weak maximum association signals (P > 1.0E-5 for GWAS cohorts)
            continue
        if file_row.loc[f'gwas_r2'] < 0.8:
            # reject cohorts showed maximal assoc with a SNP in low LD with the causal variant that we had specified
            continue
        gwas_sum_stats_list.append(pd.read_table(file_row.loc[f'gwas_sum_stats'], sep=r'\s+'))
    if len(gwas_sum_stats_list) == 0:
        print(f'warning: no suitable GWAS summary statistics file generated')
    gwas_sum_stats_df = pd.concat(gwas_sum_stats_list)
    gwas_sum_stats_df.sort_values('P', inplace=True)
    gwas_sum_stats_df.to_csv(f'gwas_simu_risk{disease_risk}_{output_suffix}.tsv', sep='\t', header=True, index=False)
    if output_dir is None:
        output_dir = f'output_{output_suffix}'
    Path(output_dir).mkdir(exist_ok=True, parents=True)
    for i in eqtl_genetic_models:
        shutil.move(f'eqtl_simu_{i}_hsq{heritability}_{output_suffix}.tsv', output_dir)
    shutil.move(f'gwas_simu_risk{disease_risk}_{output_suffix}.tsv', output_dir)
    print(f'Merge sumstat completed, duration {datetime.now() - start_time}')
    return (os.path.join(output_dir, f'gwas_simu_{output_suffix}.tsv'),
            [os.path.join(output_dir, f'eqtl_simu_{i}_{output_suffix}.tsv') for i in eqtl_genetic_models])


def get_genes(all_gene_list, output_path, gene_per_chrom=5, chrom_list=None):
    if chrom_list is None:
        chrom_list = [str(i) for i in range(1, 23)]
    start_time = datetime.now()
    print(f'Sampling gene from chromosome start at: {start_time}')
    gene_info_df = pd.read_table(all_gene_list,
                                 dtype={'start': 'Int64', 'end': 'Int64', 'position': 'Int64', 'chrom': 'string'})
    gene_info_df.drop(index=gene_info_df[~gene_info_df['chrom'].isin(chrom_list)].index,
                      inplace=True)
    gene_info_df.drop_duplicates(subset='center_snp', inplace=True)
    gene_info_df.drop(index=gene_info_df[gene_info_df['r2_0_04_cnt'] == 0].index, inplace=True)
    gene_info_df.drop(index=gene_info_df[gene_info_df['r2_04_07_cnt'] == 0].index, inplace=True)
    gene_info_df.drop(index=gene_info_df[gene_info_df['r2_07_09_cnt'] == 0].index, inplace=True)
    gene_info_df.drop(index=gene_info_df[gene_info_df['r2_09_1_cnt'] == 0].index, inplace=True)
    sample_df_list = []
    for _, grp in gene_info_df.groupby('chrom'):
        if grp.shape[0] > gene_per_chrom:
            sample_df_list.append(grp.sample(gene_per_chrom))
        else:
            sample_df_list.append(grp)
    result_df = pd.concat(sample_df_list)
    result_df['chrom'] = result_df['chrom'].astype(int)
    result_df.sort_values(['chrom', 'start'], inplace=True)
    result_df.to_csv(output_path, sep='\t', header=True, index=False)
    print(f'Sampling gene from chromosome completed, duration {datetime.now() - start_time}')


def cleanup():
    os.system('rm -f *.bed*')
    os.system('rm -f *.bim')
    os.system('rm -f *.fam')
    os.system('rm -f *.frq')
    os.system('rm -f *.legend')
    os.system('rm -f *.ld')
    os.system('rm -f *.vcf*')
    os.system('rm -f *.assoc*')
    os.system('rm -f *.nosex')
    os.system('rm -f *.par')
    os.system('rm -f *.tsv')
    os.system('rm -f *.log')
    os.system('rm -f *.txt')
    os.system('rm -f *.summary')
    os.system('rm -f *phen*')
    os.system('rm -f qtl_result*')


def run(src_vcf_dir, ld_r2_stats, genetic_map, ctrl_count=20000, case_count=20000, maf_thresh=0.1,
        sample_cnt_per_chrom=5, candicate_loci_list=None, heritability=0.2, sim_chrom_list=None, disease_risk=1.1,
        eqtl_genetic_model=2):
    start_time = datetime.now()
    print(f'Simulator start at: {start_time}')
    output_id = datetime.now().strftime("%Y%m%d%H%M%S")
    selected_gene_list_file = f'gene_list_{output_id}.tsv' if candicate_loci_list is None else candicate_loci_list
    if candicate_loci_list is None:
        full_gene_list_file = ld_r2_stats
        get_genes(full_gene_list_file, selected_gene_list_file, sample_cnt_per_chrom, sim_chrom_list)
    output_dir = f'output_{ctrl_count}_{case_count}_{output_id}'
    Path(output_dir).mkdir(exist_ok=True, parents=True)
    simulate(output_id, src_vcf_dir, selected_gene_list_file, genetic_map,
             ctrl_count, case_count, maf_thresh, heritability=heritability, disease_risk=disease_risk,
             eqtl_genetic_model=eqtl_genetic_model, output_dir=output_dir)
    print(f'Simulation and analysis completed, results in {output_dir}, total duration {datetime.now() - start_time}')
    os.system(f'mv simu*.log {output_dir} || true')
    cleanup()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--src_vcf_dir', dest='src_vcf_dir',
                        default='/Volumes/HD/biodata/colocalization-tools/raw/vcf/phased_hg38/EUR',
                        help='Directory of source vcf, vcf name must be chr{num}.vcf.gz')
    parser.add_argument('--ld_r2_stats', dest='ld_r2_stats',
                        default='/Volumes/HD/biodata/colocalization-tools/raw/genecode/ld_r2_stats.tsv',
                        help='Path of ld_r2_stats file generated by r2_calc.py')
    parser.add_argument('--genetic_map', dest='genetic_map',
                        default='/Volumes/HD/biodata/colocalization-tools/raw/genetic_map/genetic_maps.b38.tar.gz',
                        help='Path of genetic map file')
    parser.add_argument('--ctrl', dest='ctrl', default=20000, type=int, help='Number of control samples to simulate')
    parser.add_argument('--case', dest='case', default=20000, type=int, help='Number of case samples to simulate')
    parser.add_argument('--maf', dest='maf', default=0.1, type=float, help='MAF threshold to pick causal variant')
    parser.add_argument('--cnt_per_chr', dest='gene_cnt_on_chrom', default=20, type=int, help='Gene count per chrom')
    parser.add_argument('--loci_list_file', dest='loci_list_file', default=None, type=str,
                        help='Simulation gene list, --cnt_per_chr and --chrs_to_sim will be ingored if this parameter is specified')
    parser.add_argument('--heritability', dest='heritability', default=0.2, type=float,
                        help='Heritability for eqtl simulation')
    parser.add_argument('--disease_risk', dest='disease_risk', default=1.1, type=float,
                        help='Disease risk ratio for a single effect allele')
    parser.add_argument('--chrs_to_sim', dest='chroms_to_sim', default=None, type=str,
                        choices=[str(i) for i in range(1, 23)], nargs='*',
                        help='Only simulate on specified chromosomes')
    parser.add_argument('--eqtl_genetic_model', dest='eqtl_genetic_model', default=2, type=int,
                        choices=[0, 2, 3, 4],
                        help='eQTL genetic model, one of [0, 2, 3, 4], 0: h0, 2: h2004, 3:h20407, 4:h20709, '
                             'eQTL h1 (1) data will always be generated')
    args = parser.parse_args()
    print(f'Accepted args:\n {args}')
    run(args.src_vcf_dir, args.ld_r2_stats, args.genetic_map, args.ctrl, args.case, args.maf, args.gene_cnt_on_chrom,
        args.loci_list_file, args.heritability, args.chroms_to_sim, args.disease_risk, args.eqtl_genetic_model)
