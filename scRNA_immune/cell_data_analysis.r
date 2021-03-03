#Rcript
#Date: 2021-1-5
#Author: Fei Qiu



##cell data analysis
setwd("/share/pub/qiuf/COVID/01-data/CELL/dataset/E-MTAB-9357.processed.5")
file=list.files()
data.10x = list(); # first declare an empty list in which to hold the feature-barcode matrices
scrna.list = list(); 
E9357.sdrf <- read.delim("/share/pub/qiuf/COVID/01-data/CELL/dataset/E-MTAB-9357.sdrf.txt")
clinic_info_covid <- read.delim("/share/pub/qiuf/COVID/01-data/CELL/dataset/clinic_info.txt")
clinic_info_health<-read.delim("/share/pub/qiuf/COVID/01-data/CELL/dataset/clinic_info_health.txt")
infor<-data.frame(E9357.sdrf$Derived.Array.Data.File,E9357.sdrf$Factor.Value.disease.,E9357.sdrf$Factor.Value.sampling.time.point.,E9357.sdrf$Characteristics.individual.)
infor[,1]<-gsub(".gz","",infor[,1])
infor<-cbind(infor,clinic_info_covid[match(infor$E9357.sdrf.Characteristics.individual.,clinic_info_covid$Study.Subject.ID),])
clinic_info_health$Sample.ID<-paste("Healthy_",clinic_info_health$Sample.ID,sep="")
infor[0:16,-29:-1]<-clinic_info_health[!is.na(match(infor$E9357.sdrf.Characteristics.individual.,clinic_info_health$Sample.ID)),-3:-1]

print("step1: Data prepration ##################################################################")
for(idx in 1:length(file)){
  #### step1.1: create Seurat object
  data.10x[[idx]] <- read.delim( file[idx],header=TRUE,row.names = 1); 
  data.10x[[idx]]<-as(t(data.10x[[idx]]), "dgCMatrix")
  scrna.list[[idx]] = CreateSeuratObject(counts = data.10x[[idx]], min.cells=3, min.features=200, project=file[idx]);
  
  #### step1.2: mitochondrial percentage and rRNA percentage
  scrna.list[[idx]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[idx]], pattern = "^mt-")

  #### QC cutoff
  ## Seurat pipeline
  ## cutoff from [Monika Litviňuková, et al. Nature, 2020]
  scrna.list[[idx]] <- subset(scrna.list[[idx]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & nCount_RNA < 10000 )
  
  print(scrna.list[[idx]])
  
   #### step1.3: calculate cell cycle gene score
  s.genes = cc.genes.updated.2019$s.genes
  g2m.genes = cc.genes.updated.2019$g2m.genes
  scrna.list[[idx]] <- NormalizeData(scrna.list[[idx]], normalization.method = "LogNormalize", scale.factor = 10000)
  scrna.list[[idx]] <- CellCycleScoring(scrna.list[[idx]],  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  #### normalization: 'nFeature_RNA', 'nCount_RNA',"percent.mt","percent.ribo","S.Score", "G2M.Score"
  scrna.list[[idx]]=SCTransform(scrna.list[[idx]], vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
  ####get meta.data
  meta.data<-infor[match(file[idx],infor[,1]),2:length(infor[1,])]
  scrna.list[[idx]]@meta.data<-cbind(scrna.list[[idx]]@meta.data,meta.data)
}
#### step2: merge batches, add meta-information
print("step2: merge samples, add meta-information ##################################################################")
pancreas_merged <- merge(scrna.list[[1]], y = scrna.list[2:length(scrna.list)], project = "merged", merge.data = TRUE)


#pancreas_merged<-readRDS("/share/pub/qiuf/COVID/01-data/CELL/cell_pancreas_merged.rds")

## plot pre Harmony with PCA through integrationFeature
pancreas.features <- SelectIntegrationFeatures(object.list = scrna.list, nfeatures = 3000)
VariableFeatures(pancreas_merged) <- pancreas.features
pancreas_merged <- RunPCA(object = pancreas_merged, assay = "SCT", features = pancreas.features, npcs = 20)
#### step3: Harmony
pancreas_merged <- RunHarmony(object = pancreas_merged,
                              assay.use = "SCT",
                              reduction = "pca",
                              dims.use = 1:20,
                              group.by.vars = "orig.ident",
                              plot_convergence = TRUE)

#### step4: cell cluster
pancreas_merged <- RunUMAP(object = pancreas_merged, assay = "SCT", reduction = "harmony", dims = 1:20)
pancreas_merged <- FindNeighbors(object = pancreas_merged, assay = "SCT", reduction = "harmony", dims = 1:20)
pancreas_merged <- FindClusters(pancreas_merged,resolution =0.5, verbose = FALSE)
#pancreas_merged<-readRDS("/share/pub/qiuf/COVID/01-data/CELL/Cluster_0.5.rds")
options(future.globals.maxSize= 891289600)
statas<-pancreas_merged@meta.data$Who.Ordinal.Scale
statas[pancreas_merged@meta.data$Who.Ordinal.Scale%in%c("1","1 or 2")]<-"mild"
statas[pancreas_merged@meta.data$Who.Ordinal.Scale%in%c("3","4")]<-"moderate"
statas[pancreas_merged@meta.data$Who.Ordinal.Scale%in%c("5","6","7")]<-"severe"
statas[is.na(pancreas_merged@meta.data$Who.Ordinal.Scale)]<-"normal"
pancreas_merged@meta.data$statas<-statas
#### step5: cell annotation
id<-c(0:(length(unique(pancreas_merged@meta.data$seurat_clusters))-1))
marker_list<-get_marker(pancreas_merged,id)
all <- marker_list[[1]] %>% mutate(avg_fc = (severe_avg_log2FC + mild_avg_log2FC+normal_avg_log2FC+moderate_avg_log2FC) /4) %>% arrange(desc(avg_fc)) 
all<-all %>% arrange(cluster) 
all$anno<-marker[match(all$gene,marker$V2),1]
if(!is.null(marker_list[[2]])){
  all2<-marker_list[[2]] %>% group_by(cluster) %>% arrange( desc(avg_log2FC)) 
  all2<-all2 %>% arrange(cluster) 
  all2$anno<-marker[match(all2$gene,marker$V2),1]
}

annotation<-read.delim("/share/pub/qiuf/COVID/03-procesure/cell_data_processure/annotation.txt",header=FALSE)
pancreas_merged$annotation<-annotation[match(pancreas_merged$seurat_clusters,annotation$V1),3]
DimPlot(pancreas_merged, reduction = "umap", group.by = "annotation", pt.size = .1,label = T)
###step6:recluster
setwd("/share/pub/qiuf/COVID/03-procesure/cell_data_processure")
marker <- read.delim("E:/COVID/03-procesure/Gsemerged_data_processure/marker.txt", header=FALSE)
CD8Tmarker<-recluster(pancreas_merged,anno="CD8+ T",marker,0.2)
CD4Tmarker<-recluster(pancreas_merged,anno="CD4+ T",marker,0.1)
#### step7: extract Rploy predata from pancreas_merged.rds
#use SCT counts and cluster(已处理好)
CD8T<-readRDS("/share/pub/qiuf/COVID/03-procesure/cell_data_processure/CD8T.rds")
CD4T<-readRDS("/share/pub/qiuf/COVID/03-procesure/cell_data_processure/CD4T.rds")
Effector_CD8T<-subset(CD8T,subset = seurat_clusters==0)
Effector_CD8T@meta.data$RPoly<-0
Effector_CD8T@meta.data$annotation<-"Effector_CD8T"
Memory_CD8T<-subset(CD8T,subset = seurat_clusters==1)
Memory_CD8T@meta.data$RPoly<-1
Memory_CD8T@meta.data$annotation<-"Memory_CD8T"
Naive_CD8T<-subset(CD8T,subset = seurat_clusters%in%c(2,3))
Naive_CD8T@meta.data$RPoly<-2
Naive_CD8T@meta.data$annotation<-"Naive_CD8T"
NaiveCD4T<-subset(CD4T,subset = seurat_clusters==0)
NaiveCD4T@meta.data$RPoly<-3
NaiveCD4T@meta.data$annotation<-"NaiveCD4T"
MemoryCD4T<-subset(CD4T,subset = seurat_clusters%in%c(1,2))
MemoryCD4T@meta.data$RPoly<-4
MemoryCD4T@meta.data$annotation<-"MemoryCD4T"
CD14_monocyte<-subset(pancreas_merged,subset = annotation=="CD14+ monocyte")
CD14_monocyte@meta.data$RPoly<-5
NK<-subset(pancreas_merged,subset = annotation=="NK")
NK@meta.data$RPoly<-6
NaiveB<-subset(pancreas_merged,subset = annotation=="Naive B")
NaiveB@meta.data$RPoly<-7
CD16_monocyte<-subset(pancreas_merged,subset = annotation=="CD16+ monocyte")
CD16_monocyte@meta.data$RPoly<-8
DC<-subset(pancreas_merged,subset = annotation=="DC")
DC@meta.data$RPoly<-9
Platelet<-subset(pancreas_merged,subset = annotation=="Platelet")
Platelet@meta.data$RPoly<-10
CD34Progenitor<-subset(pancreas_merged,subset = annotation=="CD34+ Progenitor")
CD34Progenitor@meta.data$RPoly<-11
MatureB<-subset(pancreas_merged,subset = annotation=="Mature B")
MatureB@meta.data$RPoly<-12
Rpoly<-merge(Effector_CD8T,list(Memory_CD8T,Naive_CD8T,NaiveCD4T,MemoryCD4T,CD14_monocyte,NK,NaiveB,CD16_monocyte,DC,Platelet
                                ,CD34Progenitor,MatureB))
status<-unique(Rpoly@meta.data$statas)
for (j in 1:length(status)){
  temp<-subset(Rpoly,subset=statas==status[j])
  pancreas_exp<-temp@assays$SCT@counts
  pancreas_exp@Dimnames[[2]]<-paste0("Cluster",temp@meta.data$RPoly)
  exp<-as_matrix(pancreas_exp)
  #First standardize the expression profile,
  #and then calculate the average expression in the cluster as the expression amount of the cluster
  for(i in 0:12){
    print(i)
    y<-as.matrix(apply(exp,1,myfun<-function(x){sum(x[pancreas_exp@Dimnames[[2]]==paste0("Cluster",i)])}))
    #x<-apply(pancreas_exp,1,myfun<-function(x){sum(x[pancreas_exp@Dimnames[[2]]==paste0("Cluster",i)])})
    #z<-paste0("Cluster",i)
    if(i==0){
      Rpoly_exp<-data.frame(Cluster0=y)
    }else{
      Rpoly_exp<-cbind(Rpoly_exp,y)
    }
  }
  colnames(Rpoly_exp)<-paste0("Cluster",c(0:12))
  Rpoly_exp<-NormalizeData(Rpoly_exp)
  write.table(Rpoly_exp,paste("/share/pub/qiuf/COVID/03-procesure/cell_data_processure/Rploy_",status[j],"_cell.txt",sep=""))
}

#### step8:Plot celltype diffirentiated in disease
###picture the cell types enriched in different disease stages
all.freq = table(Rpoly@meta.data$statas, Rpoly@meta.data$annotation)
all.freq<-all.freq[c(3,1,2,4),]#Adjust the order of disease states,Normal,mild,moderate,severe
gse.prop<-as.data.frame(prop.table(all.freq))
colnames(gse.prop)<-c("status","CellType","proportion")
ggplot(gse.prop,aes(status,proportion,fill=CellType))+
  geom_bar(stat = "identity",position = "fill",width = 0.6)+
  ggtitle("")+theme_bw()+scale_fill_manual(values=c("#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#FFD700","#00008B","#B22222","#CCCC99","#FF0000","#FF6600"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),
        text=element_text(size=16,  family="serif"),axis.title.x=element_text(size=20,face="bold"),
        axis.title.y=element_text(size=20,face="bold"),
        axis.text.x=element_text(size=16,face="bold"),axis.text.y=element_text(size=16,face="bold"))+
  scale_y_continuous(expand = c(0,0))
###get pval by chisp.test,
#       monocyte   non-monocyte
#normal 100           7657
#mild   456           3259
p<-NULL
for(i in 1:(length(all.freq[,1])-1)){
  table<-all.freq[c(i,i+1),]
  l1<-sum(table[1,])
  l2<-sum(table[2,])
  temp<-t(as.matrix(apply(table,2,myfun<-function(x){chisq.test(matrix(c(x[1],x[2],l1,l2),nrow=2,ncol=2))$p.value})))
  name<-paste(row.names(table),collapse = "_")
  row.names(temp)<-name
  p<-rbind(p,temp)
}
all.freq<-as.data.frame.matrix(all.freq)
all.pct<-100*all.freq/rowSums(all.freq)
plot<-list()
for(i in 1:(length(all.freq[,1])-1)){
  temp_p<-as.matrix(p[i,])
  #colnames(temp_p)<-row.names(p)[i]
  colnames(temp_p)<-row.names(all.freq)[i]
  plot[[i]] = matrix_barplot(data=all.pct[c(i,i+1),], group_by=row.names(all.freq)[c(i,i+1)],pvals=temp_p,colors=set.colors)
}
#plot(plot[[1]])get the picture of matrix_barplot
