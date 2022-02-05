#remove.packages("CellChat")
#devtools::install_local("/Rgithub/CellChat-master")

setwd("/share/pub/huangyk/COVID-19/cellchat")
seurat_object <- readRDS("01-data/GSE158055/lung/Rpoly.rds")
#names(seurat_object@meta.data)
#head(seurat_object@meta.data$orig.ident)
#ls()
#head(seurat_object@meta.data$statas)
#seurat_object@meta.data$"percent.mt"[1:10,]
#seurat_object@meta.data$"percent.mt2"[1:10,]
#seurat_object <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt2")
library(CellChat)
library(Seurat)
doFuture::registerDoFuture()
future::plan("multiprocess", workers = 10)
options(future.globals.maxSize = 8000 * 1024^2)

run_cell_chat_seurat_obj <- function(seurat_object,output){
    data.input <- seurat_object@assays$RNA@data

    #name<-c(0:12)
    #name1<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell","Naive CD4+T cell","Memory CD4+T cell ","CD14+monocyte","NK","Naive B cell","CD16+monocyte","DC","Platelet"
    #        ,"CD34+Progenitor","Mature B cell")
    #seurat_object@meta.data$annotation<-name1[match(seurat_object@meta.data$RPoly,name)]

    labels <- seurat_object@meta.data$anno
    meta <- data.frame(group = labels, row.names = colnames(data.input))
    #name<-c(0:12)
    #name1<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell","Naive CD4+T cell","Memory CD4+T cell ","CD14+monocyte","NK","Naive B cell","CD16+monocyte","DC","Platelet"
    #         ,"CD34+Progenitor","Mature B cell")
    #meta$group <- name1[match(meta$group,name)]

    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    CellChatDB.use <- CellChatDB
    cellchat@DB <- CellChatDB.use

    cellchat <- subsetData(cellchat)

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)

    cellchat <- computeCommunProb(cellchat)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat,file=output)
}

seurat_object<- subset(x = seurat_object, subset = anno != "unknow")
name_new <- c(
        "Memory CD8+ T cell","Effector CD8+ T cell","Proliferating CD8+ T cell","Effector CD4+ T cell","Exhaustion CD4+ T cell",
        "Memory CD4+ T cell","Proliferating CD4+ T cell","B cell","Plasma B cell","CD14+ CD16+ monocyte","CD16+ monocyte","Macrophage","DC","Epithelial cell",
        "NK","Neutrophil","Plasma cell","Mast cell"
         )
name_raw <- c(
        "Memory CD8+T cell","Effector CD8+T cell","Proliferating CD8+T cell","Effector CD4+T cell","Exhaustion CD4+T cell",
        "Memory CD4+T cell" ,"Proliferating CD4+T cell","B","Plasma B","CD14+CD16+monocyte","CD16+monocyte","Macro","DC","Epi",
        "NK","Neu","Plasma","Mast"
         )
seurat_object@meta.data$anno <- name_new[match(seurat_object@meta.data$anno,name_raw)]

unique(seurat_object@meta.data$anno)
CCR1<-seurat_object@assays$SCT@counts[match("CCR1",row.names(seurat_object)),]
positive<-CCR1>0
negative<-CCR1==0
seurat_object@meta.data$anno[(seurat_object@meta.data$anno =="CD16+ monocyte")& positive] <- "CCR1+ CD16+ monocyte"
seurat_object@meta.data$anno[(seurat_object@meta.data$anno =="CD16+ monocyte")& negative] <- "CCR1- CD16+ monocyte"
seurat_object@meta.data$anno[seurat_object@meta.data$anno =="Plasma B cell"] <- "B cell"

CXCR6 <- seurat_object@assays$SCT@counts[match("CXCR6",row.names(seurat_object)),]
positive<-CXCR6>0
negative<-CXCR6==0
seurat_object@meta.data$anno[(seurat_object@meta.data$anno =="Memory CD8+ T cell")& positive] <- "CXCR6+ Memory CD8+ T cell"
seurat_object@meta.data$anno[(seurat_object@meta.data$anno =="Memory CD8+ T cell")& negative] <- "CXCR6- Memory CD8+ T cell"
unique(seurat_object@meta.data$anno)

#run_cell_chat_seurat_obj(seurat_object,"lung_sub.rds")
seurat_object_mild <- subset(x = seurat_object, subset = statas == "mild/moderate")
run_cell_chat_seurat_obj(seurat_object_mild,"lung_mild.rds")
seurat_object_severe <- subset(x = seurat_object, subset = statas == "severe/critical")
run_cell_chat_seurat_obj(seurat_object_severe,"lung_servere.rds")
