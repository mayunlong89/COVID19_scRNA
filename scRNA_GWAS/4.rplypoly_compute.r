###########
#Compute the rolypoly result
###############
load("/share/pub/dengcy/Singlecell/COVID19/data/covid_ld.RData")
load("/share/pub/dengcy/Singlecell/COVID19/data/COVID19_GWAS_autosomes_maf.RData")
COVID19_GWAS_autosomes_maf2<-COVID19_GWAS_autosomes_maf[COVID19_GWAS_autosomes_maf$maf>0.05,]
COVID19_GWAS_Hypertension2<-COVID19_GWAS_Hypertension[COVID19_GWAS_Hypertension$maf>0.05,]
COVID19_GWAS_diabetes2<-COVID19_GWAS_diabetes[COVID19_GWAS_diabetes$maf>0.05,]


colnames(covid_ld)[7]<-"R"
lapply(unique(covid_ld$CHR_A), function(i){
  a<-data.table(covid_ld[covid_ld$CHR_A == i,])
  file_name <- paste0("/share/pub/dengcy/Singlecell/COVID19/data/LD/",i,".Rds")
  saveRDS(a, file = file_name)
})

annotation<-read.delim("/share/pub/dengcy/Singlecell/COVID19/data/annotation_cell.txt",header=F)

library("rtracklayer")
library("stringr")
gtf1 <- rtracklayer::import('/share/pub/dengcy/Singlecell/COVID19/data/gencode.v32.annotation.gtf')
gtf_df <- as.data.frame(gtf1)
rm(gtf1)
gtf_df<-gtf_df[gtf_df$type=="gene",]
geneid_df <- dplyr::select(gtf_df,c(seqnames,start,end,gene_name))
#center in TSS,set a window of 100kb
geneid_df1<-geneid_df

start<-lapply(geneid_df$start, function(x){
  ifelse(x-50000 < 1,
         x<- 1,
         x<- x-50000)
  return(x)
})
geneid_df1$start<-unlist(start)
geneid_df1$end<-geneid_df$end + 50000
b<-str_replace_all(geneid_df1$seqnames,"chr","")
geneid_df1$chrom<-b
geneid_df1<-geneid_df1[,c(5,2,3,4)]
colnames(geneid_df1)[4]<-"label"
geneid_df1<-geneid_df1[!(geneid_df1$chrom %in% c("M","X","Y")),]

##save the data
setwd("/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result")
save(geneid_df1,annotation,COVID19_GWAS_autosomes_maf2,covid_ld,file = "predata.RData")

#######################
#compute rolypoly
######################
library("rolypoly")
library("data.table")
#
index<-c("normal","mild","moderate","severe")
lapply(index,function(x){
  file_n<-paste0("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_",x,"_cell.txt")
  merge_scexpr<-read.delim(file_n,sep = " ")
  colnames(merge_scexpr)<-annotation$V2
  merge_scexpr<-merge_scexpr[apply(merge_scexpr,1,sum)!=0,]
  #create the annotation files
  gene_name<-intersect(rownames(merge_scexpr),geneid_df1$label)
  geneid_df1<-geneid_df1[geneid_df1$label %in% gene_name,]
  merge_scexpr<-merge_scexpr[gene_name,]
  geneid_df1<-geneid_df1[!duplicated(geneid_df1$label),]
#############################################
  file_na<-paste0("roly_",x,"_pre.RData")
  save(geneid_df1,merge_scexpr,file=file_na)
  })

 ld_path <- "/share/pub/dengcy/Singlecell/COVID19/data/LD"

#sim_block_annotation$label<-rownames(merge_scexprc2)[1:1000]
#Rploy_remission_GSE.txt
 rolypoly_result <- rolypoly_roll(
   gwas_data = COVID19_GWAS_autosomes_maf2,
   block_annotation = geneid_df1,
   block_data = merge_scexpr,
   ld_folder =ld_path,
   bootstrap_iters = 100
  )
  save(rolypoly_result,file = "/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/rolypoly_mild_cell.RData")


