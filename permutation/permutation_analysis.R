#Rscript
#Date: 2020-12-10
#Author: Yunlong Ma
#E-mail: glb-biotech@zju.edu.cn

#Set the work directory
#setwd("F:\\Desktop\\")
setwd("C:\\Users\\Administrator\\Desktop\\permutation_analysis")
set.seed(12345)

#Part I Read data on significant genes and background genes

#Read significant genes from MAGMA gene-based association analysis
magma_sig <- read.table("MAGMA_genes.txt", header=T)
magma_gene <- magma_sig$Gene_name

#Read significant genes from S-MultiXcan across_tissues integrative genomics analysis
multixcan_sig <- read.table("S_MultiXcan.txt", header=T)
multixcan_gene <- multixcan_sig$Gene_name

#Read significant genes from S-PrediXcan Lung integrative genomics analysis
prediXcan_lung_sig <- read.table("S_PrediXcan_lung.txt", header=T)
prediXcan_lung_gene <- prediXcan_lung_sig$Gene_name

#Read significant genes from S-PrediXcan blood integrative genomics analysis
prediXcan_blood_sig <- read.table("S_PrediXcan_blood.txt", header=T)
prediXcan_blood_gene <- prediXcan_blood_sig$Gene_name


#Read background genes of S-MultiXcan integrative genomics analysis
Backgroud_multixcan<- read.table("S_MultiXcan_background.txt", header=T)
Backgroud_multixcan_gene <- Backgroud_multixcan$Gene_name

#Read background genes of S-PrediXcan Lung integrative genomics analysis
Backgroud_PrediXcan_lung <- read.table("S_PrediXcan_lung_background.txt", header=T)
Backgroud_PrediXcan_lung_gene <- Backgroud_PrediXcan_lung$Gene_name

#Read background genes of S-PrediXcan Blood integrative genomics analysis
Backgroud_PrediXcan_blood <- read.table("S_PrediXcan_blood_background.txt", header=T)
Backgroud_PrediXcan_blood_gene <- Backgroud_PrediXcan_blood$Gene_name


#Calculate the numebr of genes in each gene set
len_Sig_magma <- length(magma_gene)
len_Sig_multixcan <- length(multixcan_gene)
len_Sig_predixcan_lung <- length(prediXcan_lung_gene)
len_Sig_predixcan_blood <- length(prediXcan_blood_gene)
len_Backgroud_multixcan_gene <- length(Backgroud_multixcan_gene)
len_Backgroud_PrediXcan_lung_gene <- length(Backgroud_PrediXcan_lung_gene)
len_Backgroud_PrediXcan_blood_gene <-length(Backgroud_PrediXcan_blood_gene)


#Part II establish a function for permutation analysis

#Permutation Function
Permut_analysis <- function(x,y,z){
  
  random_selected_genes <- sample(x,y)
  
  temp <- match(random_selected_genes,z)
  
  random_overlaped_genes <- na.omit(temp)
  
  num<-length(random_overlaped_genes)
  
  return(num)
  
}


#100000 times permutation analysis for MAGMA and S-MultiXcan analysis
X1 <- Backgroud_multixcan_gene
Y1 <- len_Sig_multixcan
Z1 <- magma_gene
results_perumt1 <- replicate(100000,Permut_analysis(X1,Y1,Z1))


#100000 times permutation analysis for MAGMA and S-PrediXcan lung analysis
X2 <- Backgroud_PrediXcan_lung_gene
Y2 <- len_Sig_predixcan_lung
Z2 <- magma_gene
results_perumt2 <- replicate(100000,Permut_analysis(X2,Y2,Z2))


#100000 times permutation analysis for MAGMA and S-PrediXcan lung analysis
X3 <- Backgroud_PrediXcan_blood_gene
Y3 <- len_Sig_predixcan_blood 
Z3 <- magma_gene
results_perumt3 <- replicate(100000,Permut_analysis(X3,Y3,Z3))


#Ploting function
Fig_random <- function(x,y,z){
  temp1 <- match(y,z)
  Observed_overlaped_genes <- na.omit(temp1)
  Observed_gene_num <- length(Observed_overlaped_genes)
  xlim_value = Observed_gene_num+20
  hist(x, col="lightblue",xlab="Counts of overlapped genes",xlim = c(0,xlim_value),main=NULL)
  abline(v=Observed_gene_num,col="red",lty="longdash")
  P_value=length(x[x>Observed_gene_num])/length(x)
  x1 = Observed_gene_num+10
  freq <- table(x)
  y1= max(freq)
  text(x1,y1,paste("P=",P_value))
  
}



freq <- table(results_perumt1)
max(freq)


#Visulization for MAGMA and S-MultiXcan analysis
Fig_random(results_perumt1,multixcan_gene,magma_gene)

#Visulization for MAGMA and S-PrediXcan Lung analysis
Fig_random(results_perumt2,prediXcan_lung_gene,magma_gene)

#Visulization for MAGMA and S-PrediXcan Blood analysis
Fig_random(results_perumt3,prediXcan_blood_gene,magma_gene)


#End






