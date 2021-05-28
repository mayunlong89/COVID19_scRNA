
#permutation analysis

#2021-05-13

getwd()
setwd('C:\\Users\\MYL\\Desktop\\')
set.seed(12345)

data1 <- read.table("expressed_genes_scRNA.txt",header = TRUE)

#set a threshold of scaled expression -0.3 for highly expressed genes or not
data1$normal[which(data1$normal> -0.3)]<-1
data1$normal[which(data1$normal< -0.3)]<-0

data1$mild[which(data1$mild> -0.3)]<-1
data1$mild[which(data1$mild< -0.3)]<-0

data1$moderate[which(data1$moderate> -0.3)]<-1
data1$moderate[which(data1$moderate< -0.3)]<-0

data1$severe[which(data1$severe> -0.3)]<-1
data1$severe[which(data1$severe< -0.3)]<-0


data2<-data1


#user-defined function

permut_analysis <- function(num){
  count <- 0
  for(i in 1:num) {
    data_test <- data2[sample(nrow(data2),size=27,replace=TRUE),]
    n1 = sum(data_test$normal)
    n2 = sum(data_test$mild)
    n3 = sum(data_test$moderate)
    n4 = sum(data_test$severe)
    p1 = n2-n1
    p2 = n3-n2
    p3 = n4-n3
    if (p1>0 & p2>0 & p3 >0){
      count=count+1 
      
  }
  
  }
  return(count)
}

permut_results <- permut_analysis(10000)
Permuted_p=permut_results/10000
Permuted_p

 





