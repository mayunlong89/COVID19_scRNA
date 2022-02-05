#!/usr/bin/env bash

#MAGMA analysis: GWAS on severe COVID-19 downloaed from COVID-19 Initiative Consortium 
#File name: COVID19_HGI_B2_ALL_leave_23andme_20201020.txt.gz.
#File location: cd /share/pub/mayl/01_COVID19_GWAS_round4/
#nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &


#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/Sherlock/MAGMA
export DATA=/share/pub/mayl/01_COVID19_GWAS_round4
export OUTPUT=/share/pub/mayl/01_COVID19_GWAS_round4/MAGMA

#Formating
#cd $DATA

#generating a location file including three Columns: SNP, CHR, POS
#   gawk '{print $13, $1, $2 }' COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt > COVID19_HGI_B2_ALL_for_MAGMA_location &

#generating a --pval file including two Columns: SNP, P
#If MAGMA detects a header in the file it will look for SNP IDs and p-values in the SNP and P column respectively. 
#If no header is found it will use the first column for SNP IDs and the second column for p-values.
#   gawk '{print $13, $9}'  COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt  > COVID19_HGI_B2_ALL_for_MAGMA.results_Pval &


#MAGMA annotation:

$MAGMA_DIR/magma \
    --snp-loc  $DATA/COVID19_HGI_B2_ALL_for_MAGMA_location  \
    --annotate window=20,20 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/COVID19_HGI_B2_ALL_for_MAGMA.hg19_SNP_Gene_annotation  


#gene-based association analysi:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/COVID19_HGI_B2_ALL_for_MAGMA.results_Pval \
    N=969689 \
    --gene-annot $OUTPUT/COVID19_HGI_B2_ALL_for_MAGMA.hg19_SNP_Gene_annotation.genes.annot \
    --out $OUTPUT/COVID19_HGI_B2_ALL_for_MAGMA.hg19_SNP_Gene_Analysis_P 

#Pathway-based analysis
#MAGMA gene set-based association analysis
$MAGMA_DIR/magma \
    --gene-results $OUTPUT/COVID19_HGI_B2_ALL_for_MAGMA.hg19_SNP_Gene_Analysis_P.genes.raw \
    --set-annot $MAGMA_DIR/KEGG_for_MAGMA_annotated.txt \
    --out $OUTPUT/COVID19_HGI_B2_ALL_for_MAGMA_Gene_set_results  




