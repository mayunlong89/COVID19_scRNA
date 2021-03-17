####################
#R: Progress the input data for MAGMA
#####################
setwd("/share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID")
COVID19_covid_filtered<-read.delim("/share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/COVID19_covid_filtered.txt")
magma_Input1<-COVID19_covid_filtered[,c("rsid","all_inv_var_meta_p","all_meta_sample_N")]
magma_Input2<-COVID19_covid_filtered[,c("rsid","CHR","POS","POS")]
write.table(magma_Input2,file="/share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/magma_Input2.txt",sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file="/share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/magma_Input1.txt",sep="\t",row.names=F,quote=F)
#
colnames(COVID19_covid_filtered)[c(13,9,11)]<-c("SNP","P","N")
write.table(COVID19_covid_filtered,file="/share/pub/dengcy/Singlecell/COVID19/MAGMA/COVID/COVID19_covid_filtered.txt",sep="\t",row.names=F,quote=F)
