library(CelliD)
library(tidyverse) 
library(ggpubr)
setDimMCSlot.Seurat <- function(X, cellEmb, geneEmb, stdev = NULL, reduction.name = "mca", assay = DefaultAssay(X), ...) {
    colnames(cellEmb) <- paste0(reduction.name, "_", seq(ncol(cellEmb)))
    colnames(geneEmb) <- paste0(reduction.name, "_", seq(ncol(geneEmb)))
    DimReducObject <- CreateDimReducObject(embeddings = cellEmb, loadings = geneEmb, key = paste0(reduction.name, "_"), assay = assay)
    X@reductions[[reduction.name]] <- DimReducObject
    if (!is.null(stdev)) {
        X@reductions[[reduction.name]]@stdev <- sqrt(stdev)
    }
    X@reductions[[reduction.name]]@misc[["mca.flag"]] <- TRUE
    return(X)
}

#seurat_object2 <- readRDS("/COVID/01-data/CELL/recluster_file/Rpoly1.rds")
#head(seurat_object@meta.data)
#rm(seurat_object)
#seurat_object <- readRDS("/COVID/01-data/CELL/Cluster_0.5.rds")
#seurat_object2 <- readRDS("/COVID/01-data/GSE158055/lung/Rpoly.rds")
#meta_data1 <- seurat_object@meta.data
#meta_data1 %>% write_tsv("meta_data1.tsv")
data("HgProteinCodingGenes")
MCA <- RunMCA.dgCMatrix(
    GetAssayData(seurat_object, "data"),
    features=rownames(seurat_object@assays$SCT@data)[rownames(seurat_object@assays$SCT@data) %in% HgProteinCodingGenes]
    )
saveRDS(MCA,"MCA.rds")

geneEmb <- MCA$featuresCoordinates
cellEmb <- MCA$cellsCoordinates
stdev <- MCA$stdev
seurat_object <- setDimMCSlot.Seurat(seurat_object, cellEmb = cellEmb, geneEmb = geneEmb, stdev = stdev, reduction.name =  "mca")

#pdf("MC_plot_all.genes_cell.pdf")
#DimPlotMC(seurat_object, reduction = "mca", group.by = "anno", features = c("LZTFL1","FYCO1","XCR1","CXCR6","CCR9","SLC6A20","CCR3","CCR1","IFNAR2",
#"CCR2","CCR5","OAS3","OAS1","SPPL2C","CRHR1","STH","MAPT","UBE2D1","KANSL1","TFAM","CCRL2",
#"IBA57","IL10RB","COX10","NCR2","FOXP4","ABO","BARHL2","CCHCR1","VSTM2A","OAS2",
#"DPP9","PGLS","DNAH3"), as.text = TRUE)
#dev.off()

all_geneset = list(Top_34_genes = c("LZTFL1","FYCO1","XCR1","CXCR6","CCR9","SLC6A20","CCR3","CCR1","IFNAR2",
"CCR2","CCR5","OAS3","OAS1","SPPL2C","CRHR1","STH","MAPT","UBE2D1","KANSL1","TFAM","CCRL2",
"IBA57","IL10RB","COX10","NCR2","FOXP4","ABO","BARHL2","CCHCR1","VSTM2A","OAS2",
"DPP9","PGLS","DNAH3"),MAGMA_genes = c("LZTFL1","FYCO1","XCR1","CXCR6","CCR9","SLC6A20","CCR3","CCR1","IFNAR2",
"CCR2","CCR5","OAS3","OAS1","SPPL2C","CRHR1","STH","MAPT","UBE2D1","KANSL1","TFAM","CCRL2",
"IBA57","IL10RB","COX10","NCR2"),
S_MultiXcan=c("LZTFL1","SLC6A20","CCR9","CXCR6","XCR1","FYCO1","CCR3","CCR1","DNAH3","CCR2","ABO","IFNAR2","IL10RB","CCR5","OAS3","CCRL2"),
`S-Predixcan_blood`=c("CCR9","PGLS"),`S-PrediXcan_lung`=c("CXCR6","CCR5","FOXP4","IL10RB","FYCO1","ABO")
)


HGT_GWAS_hint <- RunCellHGT(seurat_object, pathways = all_geneset, dims = 1:50)
seurat_object@assays[["HGT_GWAS_hint"]] <- CreateAssayObject(HGT_GWAS_hint)

seurat_object2 <- readRDS("/share/pub/qiuf/COVID/01-data/CELL/recluster_file/Rpoly1.rds")
new_data_to_add <- seurat_object2@meta.data[,47:48]
rm(seurat_object2)

saveRDS(seurat_object,"COVID19_cell_HGT.rds")
seurat_object <- readRDS("COVID19_cell_HGT.rds")

all_geneset = list(Top_34_genes = c("LZTFL1","FYCO1","XCR1","CXCR6","CCR9","SLC6A20","CCR3","CCR1","IFNAR2",
"CCR2","CCR5","OAS3","OAS1","SPPL2C","CRHR1","STH","MAPT","UBE2D1","KANSL1","TFAM","CCRL2",
"IBA57","IL10RB","COX10","NCR2","FOXP4","ABO","BARHL2","CCHCR1","VSTM2A","OAS2",
"DPP9","PGLS","DNAH3"),MAGMA_genes = c("LZTFL1","FYCO1","XCR1","CXCR6","CCR9","SLC6A20","CCR3","CCR1","IFNAR2",
"CCR2","CCR5","OAS3","OAS1","SPPL2C","CRHR1","STH","MAPT","UBE2D1","KANSL1","TFAM","CCRL2",
"IBA57","IL10RB","COX10","NCR2"),
S_MultiXcan=c("LZTFL1","SLC6A20","CCR9","CXCR6","XCR1","FYCO1","CCR3","CCR1","DNAH3","CCR2","ABO","IFNAR2","IL10RB","CCR5","OAS3","CCRL2"),
`S-Predixcan_blood`=c("CCR9","PGLS"),`S-PrediXcan_lung`=c("CXCR6","CCR5","FOXP4","IL10RB","FYCO1","ABO")
)

HGT_GWAS_hint <- RunCellHGT(seurat_object, pathways = all_geneset, dims = 1:50,p.adjust = F)
seurat_object@assays[["HGT_GWAS_hint"]] <- CreateAssayObject(HGT_GWAS_hint)

meta.data.combined <- left_join(rownames_to_column(seurat_object@meta.data), rownames_to_column(new_data_to_add))
use.cells <- meta.data.combined$rowname
meta.data.combined[,1] <- NULL
rownames(meta.data.combined) <- use.cells
seurat_object@meta.data$statas <- meta.data.combined$statas
seurat_object@meta.data$annotation <- meta.data.combined$annotation

seurat_object <- subset(x = seurat_object, subset = (statas == "normal" | statas == "mild" | statas == "moderate" | statas == "severe"))
#seurat_object@assays[["HGT_GWAS_hint"]][1:5,1:5]

pdf(glue::glue("Top_34_genes_umap.pdf"))
FeaturePlot(seurat_object, reduction = "umap",features = "Top-34-genes",  cols = c("lightgrey", "#DE2D26"),raster=FALSE,min.cutoff = 1)
dev.off()
#F5F5BA
pdf(glue::glue("MAGMA_genes_umap.pdf"))
FeaturePlot(seurat_object, reduction = "umap",features = "MAGMA-genes",  cols = c("lightgrey", "#DE2D26"),raster=FALSE,min.cutoff = 1)
dev.off()

pdf(glue::glue("S_MultiXcan_genes_umap.pdf"))
FeaturePlot(seurat_object, reduction = "umap",features = "S-MultiXcan",  cols = c("lightgrey", "#DE2D26"),raster=FALSE,min.cutoff = 1)
dev.off()
#head(seurat_object@meta.data)
#nrow(seurat_object2@meta.data[,47:48])

pdf(glue::glue("Top_34_genes_umap_normal.pdf"))
FeaturePlot(seurat_object, reduction = "umap",features = "Top-34-genes",  cols = c("#F5F5BA", "#DE2D26"),cells=rownames(seurat_object@meta.data[seurat_object@meta.data$statas=="normal",]))
dev.off()

pdf(glue::glue("Top_34_genes_umap_mild.pdf"),width=12,height=6)
FeaturePlot(seurat_object, reduction = "umap",features = "Top-34-genes",  cols = c("#F5F5BA", "#DE2D26"),cells=rownames(seurat_object@meta.data[seurat_object@meta.data$statas=="mild",]))
dev.off()

pdf(glue::glue("Top_34_genes_umap_mild.pdf"),width=12,height=6)
FeaturePlot(seurat_object, reduction = "umap",features = "Top-34-genes",  cols = c("#F5F5BA", "#DE2D26"),cells=rownames(seurat_object@meta.data[seurat_object@meta.data$statas=="moderate",]))
dev.off()

pdf(glue::glue("Top_34_genes_umap_mild.pdf"),width=12,height=6)
FeaturePlot(seurat_object, reduction = "umap",features = "Top-34-genes",  cols = c("#F5F5BA", "#DE2D26"),cells=rownames(seurat_object@meta.data[seurat_object@meta.data$statas=="severe",]))
dev.off()
