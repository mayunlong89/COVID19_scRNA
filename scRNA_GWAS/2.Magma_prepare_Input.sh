##############
#Prepare the Input files for Magma
##################
#1.Get annotation for MAGMA
#Set up 50kb and down 50kb of window for each gene
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
./magma --annotate window=50,50 --snp-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/magma_Input2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/int.annotated_50kbup_50_down

#2.Get gene-level association
./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/magma_Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/int.annotated_50kbup_50_down.genes.annot \
--out /share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/int.annotated_50kbup_50_down
#3.Prepare the stat files
source activate ldsc
/share/pub/dengcy/Singlecell/COVID19/LDSC/ldsc-master/munge_sumstats.py \
--sumstats /share/pub/dengcy/Singlecell/COVID19/LDSC/COVID/COVID19_covid_filtered.txt \
--merge-alleles /share/pub/dengcy/Singlecell/COVID19/LDSC/w_hm3.snplist \
--signed-sumstat BETA,0 \
--N-col N \
--out /share/pub/dengcy/Singlecell/COVID19/LDSC/COVID/COVID
