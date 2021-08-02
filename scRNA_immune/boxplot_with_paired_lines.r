library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggforce)
library(reshape2)
library(dplyr)

#install.packages("ggforce")

data<-readRDS("E:/COVID/liver/data_0.8.rds")
gene<-read.delim("E:/COVID/liver/gene.txt",header=F)
data<-data[,data@meta.data$anno=="Cholangiocyte"]
CCR1<-data@assays$SCT@counts[match("ORMDL3",row.names(data)),]
double_positive<-CCR1>0
double_negative<-CCR1==0
positive<-data[,double_positive]
negative<-data[,double_negative]
for(i in 1:length(gene[,1])){
  score1<-NULL
  #for(z in unique(CD16@meta.data$statas)){#??ע??
  #gene=c("IFNG","IL2")
  for(j in unique(data@meta.data$orig.ident)){
    print(j)
    tryCatch({temp<-subset(positive,orig.ident==j)
             temp<-temp@assays$SCT@counts[gene[i,1],]}
             ,error=function(e){
               temp<-0
             },finally={
               score1<-rbind(score1,cbind(mean(temp),j))})

  }
score2<-NULL
for(j in unique(data@meta.data$orig.ident)){
      print(j)
      tryCatch({temp<-subset(negative,orig.ident==j)
      temp<-temp@assays$SCT@counts[gene[i,1],]}
      ,error=function(e){
        temp<-0
      },finally={
        score2<-rbind(score2,cbind(mean(temp),j))})

    }
  #}#??ע??
  score1<-cbind(score1,"ORMDL3+")
  score2<-cbind(score2,"ORMDL3-")
  temp<-rbind(score1,score2)
  temp[,1]<-as.numeric(temp[,1])
  temp[is.na(temp)]<-0
  colnames(temp)<-c("Expression","sample","ORMDL3")
  row.names(temp)<-temp[,2]
  temp<-temp[,c(1,3)]
  temp<-data.frame(temp)
  temp$ORMDL3<-factor(temp$ORMDL3,levels = c("ORMDL3-","ORMDL3+"))

  #????ͼ
  exp<-temp
  #a<-cbind(exp[exp$stage%in%c("mild","normal"),],class="a")
  #b<-cbind(exp[exp$stage%in%c("moderate","normal"),],class="b")
  a<-cbind(exp[exp$ORMDL3%in%c("ORMDL3-","ORMDL3+"),],class="a")
  exp<-a
  exp<-exp[order(exp$ORMDL3),]
  #a<-cbind(a,index=rep(seq(1:sum(a$stage=="normal")),4))
  exp<-cbind(exp,index=rep(seq(1:sum(exp$ORMDL3=="ORMDL3+")),2))
  bezier<-data.frame(matrix(ncol=4,nrow=0))
  colnames(bezier)<-c("index","ORMDL3","class","value")
  exp<-exp[,c(1,3,2,4)]
  exp$Expression<-as.numeric(exp$Expression)
  data0<-exp
  data1<-split(data0,data0$ORMDL3)
  data2<-data.frame(`ORMDL3-`=data1$`ORMDL3-`[,1],`ORMDL3+`=data1$`ORMDL3+`[,1],
                    index=data1$`ORMDL3-`[,4])
  colnames(data2)<-c(1,2,'index')
  data2$'1.3'<-data2$`1`
  data2$'1.7'<-data2$`2`
  data3<-melt(data2,id='index',variable.name = "ORMDL3")
  data3$ORMDL3<-as.numeric(as.character(data3$ORMDL3))
  data4<-arrange(data3,index,ORMDL3)
  data4$class="a"

  bezier<-rbind(bezier,data4)
  bezier$value<-as.numeric(bezier$value)
  #bezier<-cbind(melt(exp,id="index",variable.name ="stage" ),exp$class)
  my_compare<-list(c("ORMDL3+","ORMDL3-"))
  p1<-ggplot(exp, aes(x=ORMDL3, y=Expression)) + geom_boxplot(aes(fill = ORMDL3),width=0.35,alpha=0.7,position = position_dodge(0),size=0.5,outlier.size = 0)+theme_bw()+scale_fill_manual(values=c("#527239","#C31015"))+
    geom_point(aes(fill=ORMDL3),shape=21,color="black",size=2)+
    geom_bezier(data=bezier,aes(x=ORMDL3,y=value,group=index),size=0.25,color="grey20")+theme_bw()+
    stat_boxplot(geom="errorbar",width=0.15) +
    theme(plot.title = element_text(hjust = 0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"))+
    stat_compare_means(comparisons = my_compare,method="t.test",paired =T)+xlab(gene[i,1])
  save_plot(paste0("E:/COVID/liver/plot/up/curve.",gene[i,1],".pdf"),p1,base_aspect_ratio = 1)

}


a<-as.numeric(temp[temp$ORMDL3=="ORMDL3+",1])
b<-as.numeric(temp[temp$ORMDL3=="ORMDL3-",1])

#
#+facet_grid(.~class,drop=T)
