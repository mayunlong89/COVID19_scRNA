#Rscript
#Date:2020-12-12
#Author: Yunlong Ma
#E-mail: glb-biotech@zju.edu.cn



###all the top KEGG pathways in CAD:
setwd("C:\\Users\\Administrator\\Desktop\\06_KEGG_MDS_analysis_COVID19")

if(!require("proxy"))install.packages("proxy")
library(proxy)
if(!require("reshape"))install.packages("reshape")
library(reshape)

path1=read.delim("pathway-KEGG-COVID19.txt",header=FALSE)

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

Data_11 <- read.table("Jaccard_distance_final_COVID19.txt",header = FALSE)
anno <- read.table("name_KEGG_COVID19.txt",header=TRUE)  

colnames(Data_11) <- c("time","variable","value")

Data_11$time <- anno$ID[match(Data_11$time,anno$Name)] 


data_cast<-cast(Data_11,time~variable)  

data_cast$time <- anno$Name[match(data_cast$time,anno$ID)]  


#### 
data_cast_2 <-as.dist(data_cast) ### 

voles.mds=cmdscale(data_cast_2,k=9,eig=T)

sum(abs(voles.mds$eig[1:2]))/sum(abs(voles.mds$eig))

sum((voles.mds$eig[1:2])^2)/sum((voles.mds$eig)^2)

MDS_component_1 = voles.mds$points[,1]

MDS_component_2 = voles.mds$points[,2]

plot(MDS_component_1,MDS_component_2,cex=1.8,pch=21,col = "black", bg = "green")

text(MDS_component_1,MDS_component_2, row.names(voles.mds$points),cex=0.6,pos=4) ### 



#### 

dat_1 <- read.table("anno_KEGG_COVID19.txt",header=TRUE)

data_temp <- data.frame(dat_1,MDS_component_1,MDS_component_2)


#assign color to each circle
#cPal <- colorRampPalette(c('yellow','red'))
#datCol<-cPal(18)[as.numeric(cut(dat_1$Zscore,breaks = 18))]

cPal <- colorRampPalette(c("lightblue",'blue'))

datCol<-cPal(10)[as.numeric(cut(dat_1$Genesizes,breaks = 10))]


######figure plot
symbols(MDS_component_1,MDS_component_2,circle = data_temp$Zscore,inches=0.15,col = "black", bg = datCol)

text(MDS_component_1,MDS_component_2, row.names(voles.mds$points),cex=0.6,pos=4)  

