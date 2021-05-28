#Rscript
#Date:2021-05-17
#Author: Yunlong Ma
#E-mail: glb-biotech@zju.edu.cn



###all the top KEGG pathways in CAD:
#setwd("C:\\Users\\Administrator\\Desktop\\06_KEGG_MDS_analysis_COVID19")
#setwd("F:\\Desktop\\")
setwd("C:\\Users\\MYL\\Desktop\\Jaccard_distance_for_42_Pathways")


if(!require("proxy"))install.packages("proxy")
library(proxy)
if(!require("reshape"))install.packages("reshape")
library(reshape)

#path1=read.delim("enriched_pathway_42.txt",header=FALSE) #for CCR1+CD16+Monocytes
path1=read.delim("enriched_pathway_53.txt",header=FALSE) #for ABO+Megakaryocytes

path<- path1[,-1]

row.names(path)<-path1[,1]  

path2<-as.matrix(path)

number <- seq(1,length(path[,1]),1)

# 
for (i in number) {
  
  for (j in number){
    
    d1<-intersect(path2[i,],path2[j,])
    d2<-union(path2[i,],path2[j,])
    JD=length(d1)/length(d2)
    JD1=1-JD
    dnb= row.names(path2)[i]
    dna=row.names(path2)[j]
    
    DD <- c(dna,dnb,JD)
    Y =t(DD)
    
    write.table(Y, file="Jaccard_distance_final_COVID19.txt", 
                append=TRUE, row.names=FALSE,col.names = FALSE)    
    
  }
  
}

