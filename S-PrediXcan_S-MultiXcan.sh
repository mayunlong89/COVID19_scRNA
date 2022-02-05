#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export OUTPUT=/share/pub/mayl/01_COVID19_GWAS_round4/S_MultiXcan

#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#step: GWAS parsing.py

#input file: COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $DATA/03_COVID19_GWAS_round_4th_20201020v/COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map all_inv_var_meta_beta effect_size \
-output_column_map  all_inv_var_meta_p pvalue \
-output_column_map CHR chromosome \
--chromosome_format \
-output_column_map POS position \
--insert_value sample_size 969689 --insert_value n_cases 7885 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done



#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt



