###############
#Compute the top10 genes for diffirent traits based on singel cells
#################
library(tidyverse)
library("rhdf5")
library("snow")

setwd("/share/pub/dengcy/Singlecell/COVID19")
merge_scexpr1<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_Normal_GSE.txt",sep = " ")
merge_scexpr2<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_mild_GSE.txt",sep = " ")
merge_scexpr3<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_severe_GSE.txt",sep = " ")
merge_scexpr4<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_remission_GSE.txt",sep = " ")
annotation<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/annotation.txt",header=F)
gse_expr<-lapply(list(normal=merge_scexpr1,mild=merge_scexpr2,severe=merge_scexpr3,remission=merge_scexpr4),function(expr){
colnames(expr)<-annotation$V2
expr<-na.omit(expr)
expr<-expr[apply(expr,1,sum)!=0,]
return(expr)
})

merge_scexpr5<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_mild_cell.txt",sep = " ")
merge_scexpr6<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_moderate_cell.txt",sep = " ")
merge_scexpr7<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_normal_cell.txt",sep = " ")
merge_scexpr8<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_severe_cell.txt",sep = " ")
annotation<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/annotation_cell.txt",header=F)

cell_expr<-lapply(list(mild=merge_scexpr5,moderate=merge_scexpr6,normal=merge_scexpr7,severe=merge_scexpr8),function(expr){
colnames(expr)<-annotation$V2
expr<-expr[apply(expr,1,sum)!=0,]
return(expr)
})

##########################
###1.load the gene coordinate, select the 50kb window for the up and down
############################
#Filtered to remove extended MHC (chr6, 25Mb to 34Mb).
gene_coordinates <-
  read_tsv("/share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-50000<0,0,X3-50000),end=X4+50000) %>%
  select(X2,start,end,6,1) %>%
  rename(chr="X2", Gene="X6",EntrezGene="X1") %>%
  mutate(chr=paste0("chr",chr))

#########function to compute the top10 gene
top10_function <-function(exp){
#
#exp<-gse_expr[[4]]
#exp<-exp[,-3]

exp$Gene <- rownames(exp)
#Only keep genes with a unique name and tidy data.
exp <- exp %>% add_count(Gene) %>%
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as.tibble()

#############################
###2.Each cell type is scaled to the same total number of molecules.
############################
exp <- exp %>%
  group_by(column) %>%
  mutate(Expr_sum_mean=Expr*1e6/sum(Expr))
#write_tsv(exp,"single_cell_clusterExpression.txt")
##############################
###3.Specificity Calculation
#####################################
#The specifitiy is defined as the proportion of total expression performed by the cell type of interest (x/sum(x)).

exp<- exp %>%
  group_by(Gene) %>%
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>%
  ungroup()
###################
###4.Get MAGMA genes
######################
#Only keep genes that are tested in MAGMA
exp2 <- inner_join(exp,gene_coordinates,by="Gene")
##################
###5.Get number of genes
#################
#Get number of genes that represent 10% of the dataset

n_genes <- length(unique(exp2$EntrezGene))
n_genes_to_keep <- (n_genes * 0.1) %>% round()
##################
###7.Write MAGMA input files
#########################
#Filter out genes with expression below 1 TPM.
#exp3<-exp2
exp2 %>% filter(Expr_sum_mean>1) %>% magma_top10("column")
#exp3 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("column")
print("sucess!")
}

##########################
###Functions
#####################
##Get MAGMA input top10%
magma_top10 <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity)
  d_spe %>% do(write_group_magma(.,Cell_type))
}
write_group_magma  = function(df,Cell_type) {
  df <- select(df,column,EntrezGene)
  df_name <- make.names(unique(df[1]))
  colnames(df)[2] <- df_name  
  dir.create(paste0("MAGMA/"), showWarnings = FALSE)
  select(df,2) %>% t() %>% as.data.frame() %>% rownames_to_column("Cat") %>%
    write_tsv("/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10.txt",append=T)
  return(df)
}

#####run the top10_function
top10_function(gse_expr[[1]])
lapply(gse_expr,top10_function)
lapply(cell_expr,top10_function)
###spilit the top10 result files
top10_gse<-read.delim("/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_gse.txt",header=F)
gse_mild<-read.delim("/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_gse_mild.txt",header=F)
gse_normal<-read.delim("/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_gse_normal.txt",header=F)
gse_remission<-read.delim("/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_gse_remission.txt",header=F)
gse_severe<-read.delim("/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_gse_severe.txt",header=F)

top10_cell<-read.delim("/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_cell.txt",header=F)

mild=top10_cell[1:13,]
moderate=top10_cell[14:26,]
normal=top10_cell[27:39,]
severe=top10_cell[40:52,]

write.table(mild,file="/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_cell_mild.txt",quote=F)
write.table(moderate,file="/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_cell_moderate.txt",quote=F)
write.table(normal,file="/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_cell_normal.txt",quote=F)
write.table(severe,file="/share/pub/dengcy/Singlecell/COVID19/MAGMA/top10_cell_severe.txt",quote=F)
