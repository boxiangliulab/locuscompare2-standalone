# Colotools

Colotools is an interactive visualization tool for genetic association analysis of GWAS dataset and eQTL dataset.

## Overview

Colotools integrates 5 popular colocalization tools:
+ Loci-level colocalization: [JLIM](https://github.com/cotsapaslab/jlim)
+ SNP-level colocalization: [coloc](https://github.com/chr1swallace/coloc) and [fastEnloc](https://github.com/xqwen/fastenloc)
+ Mendelian randomization: [SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview)
+ Transcriptomic association: [PrediXcan](https://github.com/hakyimlab/PrediXcan) and [TWAS](http://gusevlab.org/projects/fusion/)

It could run all the colocalization tools above, display summary report and give manhattan plot and LocusCompare plot for 
the significant SNPs and genes.

## 1. Setup Environment
We provide a shell to set up the running environment for colotools.
1) Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
2) Clone colotools from [github](colotools github repo)
3) Execute the [environment set up script](./running_env/setup_env.sh) in /running_env folder
```shell
cd running_env
./setup_env.sh env.yml
```
4) The console shows 'All finished' when the setup is done. Activate the colotools virtual environment by execute:
```shell
conda activate colotools
```

## 2. Setup Configuration

### Build colotools config file ([Sample](./config.yml))
The colotools needs a config file to indicate the input data file, output file path, the GWAS and eQTL field name 
mapping etc.

Specification of config file example:
```yaml
# Required. The output root dir
working_dir: '/Volumes/HD/biodata/colotools-tools'
input:
  gwas:
    # Required. Input GWAS file path, the position should base on hg38 or hg19
    file: '/raw/Eczema/EAGLE_AD_GWAS_results_2015_hg38.tsv.gz'
    # Required. Trait name
    trait: 'eczemas'
    # Required. GWAS sample size
    sample_size: 116863
    # Required. Study type, cc or quant
    type: cc
    # Tell colotools the field name in input GWAS file.
    col_name_mapping:
      # Optional. The rs id field name in input GWAS file. If the input GWAS file does not contain rsid, remove this line
      snp: 'rsID'
      # Required. The chromosome field name in input GWAS file.
      chrom: 'chr'
      # Required. The snp position filed name in input GWAS file.
      position: 'hm_pos'
      # Required. The beta field name in input GWAS file. 
      beta: 'beta'
      # Required. The effect allele field name in input GWAS file. 
      effect_allele: 'reference_allele'
      # Required. The other allele field name in input GWAS file. 
      other_allele: 'other_allele'
      # Required. The p-value field name in input GWAS file. 
      pvalue: 'p.value'
      # Required. The se field name in input GWAS file. 
      se: 'se'
      # Required. The effect allele frequency field name in input GWAS file. 
      eaf: 'eaf'
  eqtl:
    # Required. Input eQTL file path, the position should base on hg38
    file: '/raw/eqtl/Spleen.tsv.gz'
    # Required. Tissue name
    tissue: 'Spleen'
    # Required. eQTL sample size
    sample_size: 147
    # Tell colotools the field name in the input eQTL file.
    col_name_mapping:
      # Optional. The rs id field name in input eQTL file. If the input eQTL file does not contain rsid, remove this line
      snp: 'rsid'
      # Required. The chromosome field name in input eQTL file.
      chrom: 'chromosome'
      # Required. The position field name in input eQTL file.
      position: 'position'
      # Required. The beta field name in input eQTL file.
      beta: 'beta'
      # Required. The alter allele field name in input eQTL file.
      alt: 'alt'
      # Required. The reference allele field name in input eQTL file.
      ref: 'ref'
      # Required. The p-value field name in input eQTL file.
      pvalue: 'pvalue'
      # Required. The se field name in input eQTL file.
      se: 'se'
      # Required. The gene id field name in input eQTL file.
      gene_id: 'molecular_trait_id'
      # Required. The minor allele frequency field name in input eQTL file.
      maf: 'maf'
      
  # Required if you want to run SMR/eCaviar/TWAS. The vcf files from 1000genomes.
  # NOTE that if you want to run TWAS, one more step: copy resource/convert_vcf_to_plink_binary.sh to vcf directory and run it to generate TWAS LDREF files
  # hg38 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/
  # hg19 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
  # We split the vcf files by population EUR, EAS, SAS, AMR, AFR. Please download your research population vcf files from:
  # EUR:
  # EAS:
  # SAS:
  # AMR:
  # AFR：
  # After download the vcf files, set the file directory path here
  vcf: '/Volumes/HD/biodata/colocalization-tools/raw/vcf/hg38'
  
  # Required if you want to run SMR.
  # If you want to use GTEx eQTL, download https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.basic.annotation.gtf.gz for GTEx V8,
  # or download https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz for GTEx V7.
  # Then set the file's local path here.
  # If you want to use your own eQTL, set the eQTL reference gencode version path here.
  genecode: '/Volumes/HD//biodata/colocalization-tools/raw/genecode/gencode.v26.basic.annotation.gtf.gz'
  
  # Required if you want to run fastEnloc.
  # The LD blocks are population-specific, you could download LD block file for EUR, EAS and AFR from:
  # hg38
  # EUR: https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/raw/ld_block/eur_hg38_ld_block.bed
  # EAS: https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/raw/ld_block/eas_hg38_ld_block.bed
  # AFR: https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/raw/ld_block/afr_hg38_ld_block.bed
  # hg19
  # EUR: https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/raw/ld_block/eur_hg19_ld_block.bed
  # EAS: https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/raw/ld_block/eas_hg19_ld_block.bed
  # AFR: https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/raw/ld_block/afr_hg19_ld_block.bed
  # Then set the LD block reference path here.
  ld_block_loci_file: '/Volumes/HD/biodata/colocalization-tools/raw/loci/eur_ld.hg38.bed'
  
  # Required if you want to run fastEnloc.
  # How to derive data:
  # * You could use the pre-computed GTEx V8 multi-tissue eQTL annotation here, download from https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P
  # * If you have your own eQTL data or other version of GTEx, please refer https://github.com/xqwen/fastenloc/tree/master/tutorial#derive-annotations-based-on-your-own-eqtl-data 
  #   to derive the eQTL annotation file.
  # Then set the eQTL annotation path here.
  eqtl_finemapping_file: '/Volumes/HD/biodata/colocalization-tools/raw/eqtl/finemapping/gtex_v8.eqtl_annot.vcf.gz'
  
  # Required if you want to run PrediXcan.
  # How to derive data:
  # * PrediXcan provides the prediction model and covariance file by their prediction strategies on GTEx v8 release data.
  #   You could download and unzip the data from https://zenodo.org/record/3518299/files/mashr_eqtl.tar
  # * If you have your own eQTL data or other version of GTEx, please refer https://github.com/hakyimlab/PredictDB-Tutorial to build the prediction 
  #   model and covariance.
  # Then set the prediction model and covariance based directory here.
  # Colotools will find the model and covariance in this directory by the configured eQTL tissue name.
  # The model and covariance file name format should be 'mashr_{tissue_name}.db' and 'mashr_{tissue_name}.txt.gz'.
  # For example, if the eQTL name is 'Cells_EBV-transformed_lymphocytes', Colotools will find 'mashr_Cells_EBV-transformed_lymphocytes.db' 
  # and 'mashr_Cells_EBV-transformed_lymphocytes.txt.gz' in this directory.
  prediction_dir: '/Volumes/HD/biodata/colocalization-tools/raw/prediction_model_covariance'

  # Required if you want to run JLIM.
  # How to derive data:
  # The reference genotype panels file is population-specific. Download and unpack the LD reference file by your research population: 
  # hg38:
  # EUR: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.nfe.b38.tar.gz
  # EAS: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.eas.b38.tar.gz
  # SAS: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.sas.b38.tar.gz
  # AMR: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.amr.b38.tar.gz
  # AFR：http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.afr.b38.tar.gz
  # hg19:
  # EUR: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.nfe.b37.tar.gz
  # EAS: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.eas.b37.tar.gz
  # SAS: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.sas.b37.tar.gz
  # AMR: http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.amr.b37.tar.gz
  # AFR：http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/refld.1kg.afr.b37.tar.gz
  # Then set the reference genotype panels file here.
  reference_genotype_panel_ld_ref_file: '/Volumes/HD/biodata/colocalization-tools/raw/ldref/refld.1kg.nfe.b38'
  # Required if you want to run TWAS
  # TWAS weight files of GTEx v8: http://gusevlab.org/projects/fusion/#gtex-v8-multi-tissue-expression, download and unpack the files.
  # If you want to compute your weights, refer to predictive-model-pipeline module
  twas_weights_pos: '/Volumes/HD/biodata/colocalization-tools/raw/twas_model/GTExv8.ALL.Spleen.nofilter.pos'

p-value_threshold:
  # GWAS and eQTL significance P-value threshold. Coefficient type should be float. 
  # GWAS P-value threshold must <=1.0E-7, eQTL P-value threshold must <=1.0E-5, values out of this range will be discarded
  # For example, use 5.0E-8 rather than 5E-8
  gwas: 5.0E-8
  eqtl: 1.0E-6
  
# The upstream and downstream variants number of a significant snp
neighbour_snp_range: 25
# Your research population, EUR, EAS, SAS, AMR or AFR
population: 'EUR'
```

### Generate Multiple Configs GWAS and eQTL (Optional)
If you need to run a number of GWAS and eQTL data, you need to provide config file for each GWAS-eQTL pair.

The [config generator](/common/config_generator.py) is a tool to generate config files.

+ Use Config Generator:

1. Create a GWAS config list yaml file. [Sample](/resource/gwas_config.yml)
2. Create a eQTL config list yaml file. [Sample](/resource/eqtl_config.yml)
3. Create a global config yaml file. [Sample](/resource/config_template.yml)
4. Run python script
```python
python3 colotools_project_path/common/config_generator.py --out output_dir 
[--global_cfg global_config_path] [--gwas_cfg gwas_configs_path] [--eqtl_cfg eqtl_configs_path]
```
5. Config files for each GWAS-eQTL pair will be created to the output directory.

### Config Arguments for Each Tool (Optional)

You can set arguments for each tools by creating a config file named 'tools_config.yml' in project root folder 'colocalization-tools/'.

The configuration example is in [/resource/tools_config.yml](/resource/tools_config.yml) .

## 3. Run Colotools

+ Run Colotools in command line
```shell
python3 colotools --config config_yml_file [--tools tools_names] 
```
| Parameter | Description                                                                                                                                                                                                                                     |
|-----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| config    | Required. The config file path.                                                                                                                                                                                                                 |
| tools     | Optional. The tools you want to run, support run multiple tools just split the tools name with space (e.g. --tools coloc srm jlim). Tools name could be **coloc**, **fastenloc**, **smr**, **predixcan**, **jlim**. Run all of them by default. |

+ Run Colotools in python project

```python
import colotools

colotools.run(config_yml_file, ['coloc','fastenloc']) # The tools name list. Run all tools by default.
```

### Output

After running Colotools, a summarized report and plot will be generated in the working directory.

<span style="color:red">**TODO Specification of the report and plot**</span>

+ Report
  + TODO path
  + Content and description

+ Plot
  + TODO path
  + How to use

## License
<span style="color:red">**TODO Choose an appropriate license.**</span>