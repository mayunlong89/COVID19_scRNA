#Rcript
#Date: 2020-12-10
#Author: Fei Qiu


recluster<-function(data,anno,marker,resolution){
  #data is  singlecell dataset
  #anno is the name of big-cluster
  #marker is collected from three paper
  temp<-subset(data,subset= annotation==anno)
  temp <- RunUMAP(object = temp, assay = "SCT", reduction = "harmony", dims = 1:50)
  temp <- FindNeighbors(object = temp, assay = "SCT", reduction = "harmony", dims = 1:50)
  temp <- FindClusters(object = temp, resolution = resolution)
  #plot recluster Dimplot
  png(paste(anno,resolution,".png",sep=""),width=400,height=400,res=72*2)
  DimPlot(temp,label=TRUE)
  dev.off()
  id<-c(0:(length(unique(temp@meta.data$seurat_clusters))-1))
  marker_list<-get_marker(temp,id)
  #marker_list[[1]] <- marker_list[[1]] %>% mutate(avg_fc = (severe_avg_log2FC + mild_avg_log2FC+normal_avg_log2FC+moderate_avg_log2FC) /4) %>% arrange(desc(avg_fc)) 
 #改写severe-avg-log
  marker_list[[1]] <- marker_list[[1]] %>% mutate(avg_fc=rowMeans(marker_list[[1]][,grep("_avg_log2FC$",colnames(marker_list[[1]]))]))%>% arrange(desc(avg_fc)) 
  marker_list[[1]]<-marker_list[[1]] %>% arrange(cluster)
  marker_list[[1]]$anno<-marker[match(marker_list[[1]]$gene,marker$V2),1]
  if(!is.null(marker_list[[2]])){
    marker_list[[2]]<-marker_list[[2]] %>% group_by(cluster) %>% arrange( desc(avg_log2FC)) 
    marker_list[[2]]<-marker_list[[2]] %>% arrange(cluster) 
    marker_list[[2]]<-marker[match(marker_list[[2]]$gene,marker$V2),1]
  }
  
  saveRDS(temp,paste(anno[1],".rds",sep=""))
  return(marker_list)
}
get_marker<-function(data,id){
  j=2
  marker_list=data.frame()
  all_markers<-NULL
  for(clusterID in id){
    print(clusterID)
    tryCatch({
      conserved_markers  <- FindConservedMarkers(data, ident.1 = clusterID, grouping.var = "statas", logfc.threshold = 0.25, only.pos = TRUE, verbose = FALSE) 
    }, warning = function(w){
    }, error = function(e){
      if(j!=1){
        all_markers <- FindAllMarkers(data, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        j=1
        # group exist less than 3 cells,use FindAllMarkers replace it,and only perform one time
      }
      
    },finally = {
      if(!is.na(conserved_markers[1,])){
        conserved_markers$cluster=clusterID
        conserved_markers$gene=rownames(conserved_markers)
        marker_list = rbind.fill(marker_list, conserved_markers)
      }
    })
    
  }
  marker_list[is.na(marker_list)]<-0
  list<-list(marker_list,all_markers)
  return(list)
}

matrix_barplot = function(data, group_by=NULL, pvals=NULL, xlab='', ylab='Frequency', legend.title='Groups', colors='Paired', pos='dodge', border=NA,
                          out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, sig_only=F, do.facet=F){
  
  # Plot barplot of [M x N] matrix
  # x-axis = matrix columns (e.g. cell types)
  # y-axis = matrix values (e.g. frequencies)
  # fill = matrix rows (e.g. samples) or groups (e.g. conditions)
  
  # Arguments:
  # group.by = the group of each row
  # pvals = [G x N] matrix of p-values for each group and column
  
  # Groups (default = rows)
  if(is.null(group_by)){group_by = rownames(data)}
  if(nlevels(group_by) == 0){group_by = as.factor(group_by)}
  
  # Select significant comparisons
  if(sig_only == TRUE){
    j = apply(pvals, 2, min) <= .05
    if(sum(j) == 0){return(NULL)}
    data = data[,j,drop=F]
    pvals = pvals[,j,drop=F]
  }
  
  # Construct input data
  names = colnames(data)
  data = data.frame(group=group_by, data)
  group_levels = levels(group_by)
  colnames(data)[2:ncol(data)] = names
  data = as.data.table(gather_(data, 'x', 'y', setdiff(colnames(data), 'group')))
  
  # Add p-values 1
  if(!is.null(pvals)){
    pvals = as.data.frame(pvals) %>% rownames_to_column('x') %>% gather(group, pval, -x) %>% as.data.table()
    setkeyv(data, c('x', 'group'))
    setkeyv(pvals, c('x', 'group'))
    data = merge(data, pvals, all=T)
    data$lab1 = ifelse(data$pval <= .001, '**', ifelse(data$pval <= .05, '*', ''))
  }
  
  if(coord_flip == TRUE){names = rev(names); group_levels=rev(group_levels)}
  data$x = factor(data$x, levels=names)    
  data$group = factor(data$group, levels=group_levels)
  
  # Get colors
  if(length(colors) == 1){colors = set.colors[1:length(group_levels)]}
  #plot data
  pos = position_dodge(0.9)
  p = ggplot(data) + geom_bar(aes(x=x, y=y, fill=group), colour=border, size=.7, stat='identity', position=pos)+ 
    scale_fill_manual(values=colors, name=legend.title) + xlab(xlab) + ylab(ylab) +
    scale_color_manual('', values=c('#000000', '#999999', '#cccccc'), guide='none')+theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.border = element_blank(),axis.line = element_line(colour = "black"),
          text=element_text(size=16,  family="serif"),axis.title.x=element_text(size=20,face="bold"),
          axis.title.y=element_text(size=20,face="bold"),
          #axis.text.x=element_text(size=13,face="bold",angle=-15,hjust=0.05),
          axis.text.y=element_text(size=13,face="bold"),
          legend.text=element_text(size=16))+
    geom_text(aes(x=x, y=y, label=lab1, group=group), vjust=1.0, size=6, angle=90, position=position_dodge(0.9))+
    scale_y_continuous(limits=c(0, 45),expand = c(0,0))+coord_flip()
  
  
  # Save plot
  #if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
  #p
}
as_matrix <- function(mat){
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
