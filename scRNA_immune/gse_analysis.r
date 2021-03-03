#Rcript
#Date: 2020-12-10
#Author: Fei Qiu



###gse dataset analysis
path1<-"E:/COVID/01-data/GSE149689"
path2<-"E:/COVID/01-data/GSE150861"
matrix1<-Read10X(path1)
matrix2<-Read10X(path2)
filter_10x_object1<-CreateSeuratObject(counts=matrix1,min.cells = 3,min.genes=200,project="GSE149689")
filter_10x_object2<-CreateSeuratObject(counts=matrix2,min.cells = 3,min.genes=200,project="GSE150861")
###step1:预处理数据
#meta标签增加,处理GSE149689
label<-unlist(filter_10x_object1@assays$RNA@counts@Dimnames[2])
label<-gsub(".{16}-","",label)
label<-gsub("^1$","severe",label)
label<-gsub("^2$","mild",label)
label<-gsub("^3$","flu",label)
label<-gsub("^4$","flu",label)
label<-gsub("^5$","Normal",label)
label<-gsub("^6$","flu",label)
label<-gsub("^7$","flu",label)
label<-gsub("^8$","flu",label)
label<-gsub("^9$","severe",label)
label<-gsub("10","severe",label)
label<-gsub("11","mild",label)
label<-gsub("12","mild",label)
label<-gsub("13","Normal",label)
label<-gsub("14","Normal",label)
label<-gsub("15","severe",label)
label<-gsub("16","severe",label)
label<-gsub("17","severe",label)
label<-gsub("18","mild",label)
label<-gsub("19","Normal",label)
label<-gsub("20","mild",label)
filter_10x_object1@meta.data$statas<-label
#提取不为flu的子集
filter_10x_object1<-subset(filter_10x_object1, subset = statas!="flu" )
#处理GSE150861
label<-unlist(filter_10x_object2@assays$RNA@counts@Dimnames[2])
label<-gsub(".{16}-","",label)
label<-gsub("[3,4,6,7]","remission",label)
label<-gsub("[1,2,5]","severe",label)
filter_10x_object2@meta.data$statas<-label
options(future.globals.maxSize = 10000 * 1024^2)
scrna.list<-list()
scrna.list[[1]]<-filter_10x_object1
scrna.list[[2]]<-filter_10x_object2

#### create Seurat object
#### mitochondrial percentage and rRNA percentage
scrna.list[[1]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[1]], pattern = "^MT-")
scrna.list[[2]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[2]], pattern = "^MT-")
#scrna.list[[1]][["percent.ribo"]] <- PercentageFeatureSet(scrna.list[[1]], pattern = "^Rp[sl][[:digit:]]")
#scrna.list[[2]][["percent.ribo"]] <- PercentageFeatureSet(scrna.list[[2]], pattern = "^Rp[sl][[:digit:]]")
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
#### QC cutoff
###step2:质控
## cutoff from [GSE paper]
scrna.list[[1]]<- subset(scrna.list[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15 & nCount_RNA > 500 )
scrna.list[[2]] <- subset(scrna.list[[2]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 )
scrna.list[[1]] <- NormalizeData(scrna.list[[1]], normalization.method = "LogNormalize", scale.factor = 10000)
scrna.list[[1]] <- CellCycleScoring(scrna.list[[1]],  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scrna.list[[2]] <- NormalizeData(scrna.list[[2]], normalization.method = "LogNormalize", scale.factor = 10000)
scrna.list[[2]] <- CellCycleScoring(scrna.list[[2]],  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scrna.list[[1]]=SCTransform(scrna.list[[1]], vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
scrna.list[[2]]=SCTransform(scrna.list[[2]], vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
#pancreas_merged<-scrna.list[[1]]
pancreas_merged <- merge(scrna.list[[1]], y = scrna.list[2:length(scrna.list)], project = "merged", merge.data = TRUE)
## plot pre Harmony with PCA through integrationFeature
pancreas.features <- SelectIntegrationFeatures(object.list = scrna.list, nfeatures = 3000)
VariableFeatures(pancreas_merged) <- pancreas.features
pancreas_merged <- RunPCA(object = pancreas_merged, assay = "SCT", features = pancreas.features, npcs = 50)
#### step3: Harmony
pancreas_merged <- RunHarmony(object = pancreas_merged,
                              assay.use = "SCT",
                              reduction = "pca",
                              dims.use = 1:50,
                              group.by.vars = "orig.ident",
                              plot_convergence = TRUE)

#### step4: cell cluster
pancreas_merged <- RunUMAP(object = pancreas_merged, assay = "SCT", reduction = "harmony", dims = 1:50)
pancreas_merged <- FindNeighbors(object = pancreas_merged, assay = "SCT", reduction = "harmony", dims = 1:50)
pancreas_merged <- FindClusters(object = pancreas_merged, resolution = clusterResolution)

#### step5: find cluster biomarkers
id<-c(0:(length(unique(pancreas_merged@meta.data$seurat_clusters))-1))
marker_list<-get_marker(pancreas_merged,id)

#### step6: annotate cluster 
all <- marker_list[[1]] %>% mutate(avg_fc = (severe_avg_log2FC + mild_avg_log2FC+Normal_avg_log2FC+remission_avg_log2FC) /4) %>% arrange(desc(avg_fc)) 
all<-all %>% arrange(cluster) 
all$anno<-marker[match(all$gene,marker$V2),1]
if(!is.null(marker_list[[2]])){
  all2<-marker_list[[2]] %>% group_by(cluster) %>% arrange( desc(avg_log2FC)) 
  all2<-all2 %>% arrange(cluster) 
  all2$anno<-marker[match(all2$gene,marker$V2),1]
}
marker <- read.delim("E:/COVID/03-procesure/Gsemerged_data_processure/marker.txt", header=FALSE)

#注释信息
#pancreas_merged<-readRDS("E:/COVID/03-procesure/Gsemerged_data_processure/gse_merged.rds")
annotation<-read.delim("E:/COVID/03-procesure/Gsemerged_data_processure/annotation.txt",header=FALSE)
pancreas_merged$annotation<-annotation[match(pancreas_merged$seurat_clusters,annotation$V1),3]
DimPlot(pancreas_merged, reduction = "umap", group.by = "annotation", pt.size = .1,label = T)
DimPlot(pancreas_merged, reduction = "umap", pt.size = .1,label = T)
### 画图
DimPlot(pancreas_merged,label=TRUE)
##################需要调整
FeaturePlot(pancreas_merged, features = c("CD8A","CD3D","IL7R","CCR7","CD3E"))
FeaturePlot(pancreas_merged, features = c("CD79A","MS4A1","CD19"))
FeaturePlot(pancreas_merged, features = c("FCGR3A","NCAM1","NKG7","GZMB"))
FeaturePlot(pancreas_merged, features = c("FCGR3A","NCAM1","NKG7","GZMB"))
FeaturePlot(temp, features = c("PPBP"))
FeaturePlot(temp, features = c("CD4","CLEC4A","CD14","FCGR3A","CD68","S100A12"))
FeaturePlot(pancreas_merged, features = c("HBB"))
###重聚类，从CD8T重复到monocyte
#CD8+T 分为记忆，增殖，耗竭，效应，毒性，初始
###clusteration需要调整
marker <- read.delim("E:/COVID/03-procesure/Gsemerged_data_processure/marker.txt", header=FALSE)
CD8Tmarker<-recluster(pancreas_merged,anno="CD8+ T",marker,0.2)
Bmarker<-recluster(pancreas_merged,anno="CD8+ T",marker,0.2)
CD4Tmarker<-recluster(pancreas_merged,anno="CD8+ T",marker,0.2)
MDmarker<-recluster(pancreas_merged,anno="CD8+ T",marker,0.2)
monocytemarker<-recluster(pancreas_merged,anno="CD8+ T",marker,0.2)

#### step7: extract Rploy predata from pancreas_merged.rds
#use SCT counts and cluster(已处理好)
CD8T<-readRDS("E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/CD8+ T.rds")
Effector_CD8T<-subset(CD8T,subset = seurat_clusters%in%c(0,2,4,5,6))
Effector_CD8T@meta.data$RPoly<-0
Memory_CD8T<-subset(CD8T,subset = seurat_clusters%in%c(1,3))
Memory_CD8T@meta.data$RPoly<-1
PPBP_PF4_CD8T<-subset(CD8T,subset = seurat_clusters==7)
PPBP_PF4_CD8T@meta.data$RPoly<-2
CD4T<-readRDS("E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/CD4+ T.rds")
NaiveCD4T<-subset(CD4T,subset = seurat_clusters%in%c(0,1,3,4,5,6))
NaiveCD4T@meta.data$RPoly<-3
EffectorCD4T<-subset(CD4T,subset = seurat_clusters==2)
EffectorCD4T@meta.data$RPoly<-4
B<-readRDS("E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/B.rds")
lg_B<-subset(B,subset = seurat_clusters%in%c(1,2,3))
lg_B@meta.data$RPoly<-5
NaiveB<-subset(B,subset = seurat_clusters==0)
NaiveB@meta.data$RPoly<-6
monocyte<-readRDS("E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/Monocyte.rds")
CD14_monocyte<-subset(monocyte,subset = seurat_clusters%in%c(0,1,2,3,5,6))
CD14_monocyte@meta.data$RPoly<-7
CD16_monocyte<-subset(monocyte,subset = seurat_clusters==4)
CD16_monocyte@meta.data$RPoly<-8
NK<-subset(pancreas_merged,subset = annotation=="NK")
NK@meta.data$RPoly<-9
Platelet<-subset(pancreas_merged,subset = annotation=="Platelet")
Platelet@meta.data$RPoly<-10
Rpoly<-merge(Effector_CD8T,Memory_CD8T,PPBP_PF4_CD8T,NaiveCD4T,EffectorCD4T,
         lg_B,NaiveB,CD14_monocyte,CD16_monocyte,NK,Platelet)
#saveRDS(Rpoly,"E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/Rpoly.rds")
#Rpoly<-readRDS("E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/Rpoly.rds")
status<-unique(Rpoly@meta.data$statas)
for (j in 1:length(status)){
  temp<-subset(Rpoly,subset=statas==status[j])
  pancreas_exp<-temp@assays$SCT@counts
  pancreas_exp@Dimnames[[2]]<-paste0("Cluster",temp@meta.data$RPoly)
  exp<-as_matrix(pancreas_exp)
  #First standardize the expression profile,
  #and then calculate the average expression in the cluster as the expression amount of the cluster
  for(i in 0:10){
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
  colnames(Rpoly_exp)<-paste0("Cluster",c(0:10))
  Rpoly_exp<-NormalizeData(Rpoly_exp)
  write.table(Rpoly_exp,paste("E:/COVID/01-data/CELL/recluster_file/Rploy_",status[j],"_GSE.txt",sep=""))
}

#### step8:Plot celltype diffirentiated in disease
###picture the cell types enriched in different disease stages
all.freq = table(Rpoly@meta.data$statas, Rpoly@meta.data$annotation)
all.freq<-all.freq[c(2,1,3,4),]#Adjust the order of disease states
gse.prop<-as.data.frame(prop.table(all.freq))
colnames(gse.prop)<-c("status","CellType","proportion")
ggplot(gse.prop,aes(status,proportion,fill=CellType))+
  geom_bar(stat = "identity",position = "fill",width = 0.6)+
  ggtitle("")+theme_bw()+scale_fill_manual(values=c("#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#FFD700","#00008B","#B22222","#CCCC99"))+
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


