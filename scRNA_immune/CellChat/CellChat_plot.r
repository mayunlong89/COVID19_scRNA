library(CellChat)
library(dplyr)

plotDataPrepare <- function(cell_chat_obj,sources.use,targets.use){
  cells.level <- levels(cell_chat_obj@idents)
  sources.use <- cells.level[sources.use]
  targets.use <- cells.level[targets.use]
  df.net <- subsetCommunication(cell_chat_obj, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  thresh = 0.05)
  df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
  source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
  source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
  if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
  }

  df.net$pval[df.net$pval > 0.05] = 0.2
  df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 0.9
  df.net$pval[df.net$pval <= 0.01] = 1.8
  df.net$prob[df.net$prob == 0] <- NA
  df.net$prob.original <- df.net$prob
  df.net$prob <- -1/log(df.net$prob)

  idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 0)
  if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T)*1.1, max(df.net$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), position)]
  }
      # rownames(df.net) <- df.net$interaction_name_2

  df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
  df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
  group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")

  df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
  df.net <- with(df.net, df.net[order(interaction_name_2),])
  df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
  cells.order <- group.names
  df.net$source.target <- factor(df.net$source.target, levels = cells.order)
  df <- df.net

  min.cutoff <- quantile(df$prob, 0,na.rm= T)
  max.cutoff <- quantile(df$prob, 1,na.rm= T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff

  df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),unique(df$source.target)))

  df <- df %>% filter(!is.na(source)) %>% filter(!is.na(target))
  df
}

filter_interaction_one_cell_type <- function(cell_chat_obj,cell_type_index,exclude_index){
    celltype_levels <- levels(cell_chat_obj@idents)
    df_source <- plotDataPrepare(cell_chat_obj,cell_type_index,c(1:length(celltype_levels))[c(-cell_type_index,-exclude_index)])
    nselect <- length(unique(df_source$source))
    print(glue::glue("{nselect} cell types selected"))
    df_target <- plotDataPrepare(cell_chat_obj,c(1:length(celltype_levels))[c(-cell_type_index,-exclude_index)],cell_type_index)
    nselect <- length(unique(df_target$target))
    print(glue::glue("{nselect} cell types selected"))
    df <- bind_rows(df_source,df_target)
    df$source.target2  <- paste(df$source, df$target, sep = " - ")

    df <- df %>% mutate(xname = ifelse(source == celltype_levels[cell_type_index],as.character(target),as.character(source)))
    df <- df %>% distinct(xname,interaction_name_2,.keep_all = TRUE)
    return(df)
}

interaction_simplified <- function(df,cell_type_index,exclude_index){
    celltype_levels <- unique(pull(df,source))
    df$source.target2  <- paste(df$source, df$target, sep = " - ")     
    df <- df %>% mutate(xname = ifelse(source == celltype_levels[cell_type_index],as.character(target),as.character(source)))
    df <- df %>% distinct(xname,interaction_name_2,.keep_all = TRUE)
    return(df)
}


PBMC_normal <- readRDS("PBMC/RDS/PBMC_normal_sub.rds")
PBMC_mild <- readRDS("PBMC/RDS/PBMC_mild_sub.rds")
PBMC_moderate <- readRDS("PBMC/RDS/PBMC_moderate_sub.rds")
PBMC_severe <- readRDS("PBMC/RDS/PBMC_severe_sub.rds")

levels(PBMC_normal@idents)

df1 <- filter_interaction_one_cell_type(PBMC_normal,2,c(1,9)) %>% mutate(conditions="normal")
df2 <- filter_interaction_one_cell_type(PBMC_mild,2,c(1,9)) %>% mutate(conditions="mild")
df3 <- filter_interaction_one_cell_type(PBMC_moderate,2,c(1,9)) %>% mutate(conditions="moderate")
df4 <- filter_interaction_one_cell_type(PBMC_severe,2,c(1,9)) %>% mutate(conditions="severe")

df<- bind_rows(df1,df2,df3,df4)
#select_row <- dplyr::setdiff(df4 %>% select(xname,interaction_name_2),df1 %>% select(xname,interaction_name_2)) %>% pull(interaction_name_2)
#select_row <- unique(as.character(select_row))

select_row <- dplyr::setdiff(df4 %>% select(xname,interaction_name_2),df1 %>% select(xname,interaction_name_2)) %>% pull(interaction_name_2)
df <- df %>% filter(interaction_name_2 %in% unique(as.character(select_row)))

color.use <- RColorBrewer::brewer.pal(n = 6, name = "Reds")
pdf("CCR1_PBMC_cond4_2.pdf",width=5,height=2.25)
angle <- 45
line.size <- 0.1
hjust.x <- 1
#+  facet_wrap(~target)
g <- ggplot(df, aes(x = interaction_name_2, y = xname,size = pval,fill = prob)) +
    geom_point(pch = 21,stroke = 0.1)  + facet_grid(~factor(conditions,levels=c("normal","mild","moderate","severe"))) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust= hjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(), strip.placement = "outside",strip.text = element_text(colour = 'black')) +
    scale_x_discrete(position = "bottom")
values <- c(0.2,0.9,1.8); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob.",barheight=5))
} else {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob."),barheight=5)
}
g <- g + theme(text = element_text(size = 5),plot.title = element_text(size=5)) +
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 3))

g
dev.off()

df1 <- filter_interaction_one_cell_type(PBMC_normal,6,c(5,9)) %>% mutate(conditions="normal")
df2 <- filter_interaction_one_cell_type(PBMC_mild,6,c(5,9)) %>% mutate(conditions="mild")
df3 <- filter_interaction_one_cell_type(PBMC_moderate,6,c(5,9)) %>% mutate(conditions="moderate")
df4 <- filter_interaction_one_cell_type(PBMC_severe,6,c(5,9)) %>% mutate(conditions="severe")
df<- bind_rows(df1,df2,df3,df4)
select_row <- dplyr::setdiff(df4 %>% select(xname,interaction_name_2),df1 %>% select(xname,interaction_name_2)) %>% pull(interaction_name_2)
df <- df %>% filter(interaction_name_2 %in% unique(as.character(select_row)))
pdf("CXCR6_PBMC_cond4_2.pdf",width=6,height=2.5)
angle <- 45
line.size <- 0.1
hjust.x <- 1
#+  facet_wrap(~target)
g <- ggplot(df, aes(x = interaction_name_2, y = xname,size = pval,fill = prob)) +
    geom_point(pch = 21,stroke = 0.1)  + facet_grid(~factor(conditions,levels=c("normal","mild","moderate","severe"))) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust= hjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(), strip.placement = "outside",strip.text = element_text(colour = 'black')) +
    scale_x_discrete(position = "bottom")
values <- c(0.2,0.9,1.8); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob.",barheight=5))
} else {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob."),barheight=5)
}
g <- g + theme(text = element_text(size = 5),plot.title = element_text(size=5)) +
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 3))

g
dev.off()

df <- df %>% filter(conditions %in% c("normal","severe"))
pdf("CXCR6_PBMC_fix.pdf",width=3.5,height=1.75)
angle <- 45
line.size <- 0.1
hjust.x <- 1
#+  facet_wrap(~target)
g <- ggplot(df, aes(x = interaction_name_2, y = xname,size = pval,fill = prob)) +
    geom_point(pch = 21,stroke = 0.1)  + facet_grid(~factor(conditions,levels=c("normal","mild","moderate","severe"))) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust= hjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(), strip.placement = "outside",strip.text = element_text(colour = 'black')) +
    scale_x_discrete(position = "bottom")
values <- c(0.2,0.9,1.8); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob.",barheight=5))
} else {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob."),barheight=5)
}
g <- g + theme(text = element_text(size = 5),plot.title = element_text(size=5)) +
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 3))

g
dev.off()


unique(pull(df2,source))
#rev
df1 <- interaction_simplified(df1,6,c(5,9)) %>% mutate(conditions="normal")
df2 <- interaction_simplified(df2,6,c(5,9)) %>% mutate(conditions="mild")
df3 <- interaction_simplified(df3,6,c(5,9)) %>% mutate(conditions="moderate")
df4 <- interaction_simplified(df4,6,c(5,9)) %>% mutate(conditions="severe")
df<- bind_rows(df1,df2,df3,df4)
select_row <- dplyr::setdiff(df4 %>% select(xname,interaction_name_2),df1 %>% select(xname,interaction_name_2)) %>% pull(interaction_name_2)
df <- df[df$interaction_name_2 %in% unique(as.character(select_row)),]

color.use <- RColorBrewer::brewer.pal(n = 6, name = "Reds")
pdf("CXCR6_4statas.pdf",width=7.5,height=2)
angle <- 45
line.size <- 0.1
hjust.x <- 1
#+  facet_wrap(~target)
g <- ggplot(df, aes(x = interaction_name_2, y =xname ,size = pval,fill = prob)) +
    geom_point(pch = 21,stroke = 0.1)  + facet_grid(~factor(conditions,levels=c("normal","mild","moderate","severe"))) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust= hjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(), strip.placement = "outside",strip.text = element_text(colour = 'black')) +
    scale_x_discrete(position = "bottom")
values <- c(0.2,0.9,1.8); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob.",barheight=5))
} else {
  g <- g + scale_fill_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
    guides(fill = guide_colourbar(barwidth = 0.2, title = "Commun. Prob."),barheight=5)
}
g <- g + theme(text = element_text(size = 5),plot.title = element_text(size=5)) +
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 3))

g
dev.off()

devtools::install_local("/share/pub/huangyk/COVID-19/cellchat/CellChat-master",dependencies = F)
mamba install -c conda-forge r-svglite
mamba install -c conda-forge r-rsvg

###select_row prepare
df1_source <- plotDataPrepare(PBMC_normal,6,c(1:5,8,10:15))
df1_target <- plotDataPrepare(PBMC_normal,c(1:5,8,10:15),6)
df1 <- bind_rows(df1_source,df1_target)
df1$source.target2  <- paste(df1$source, df1$target, sep = " - ") 
df1 <- df1 %>% mutate(conditions="normal")
df1 <- df1 %>% mutate(xname = ifelse(source == "CXCR6+ Memory CD8+ T cell",as.character(target),as.character(source)))
df1 <- df1 %>% distinct(xname,interaction_name_2,.keep_all = TRUE)

df2_source <- plotDataPrepare(PBMC_severe,2,c(3:8,10:15))
df2_target <- plotDataPrepare(PBMC_severe,c(3:8,10:15),2)
df2 <- bind_rows(df2_source,df2_target)
df2$source.target2 <- paste(df2$source, df2$target, sep = " - ") 
df2 <- df2 %>% mutate(conditions="severe")
df2 <- df2 %>% mutate(xname = ifelse(source == "CXCR6+ Memory CD8+ T cell",as.character(target),as.character(source)))
df2 <- df2 %>% distinct(xname,interaction_name_2,.keep_all = TRUE)
df<- bind_rows(df1,df2)
dfp <- df %>% filter(conditions=="severe")
dfn <- df %>% filter(conditions=="normal")
select_row2 <- dplyr::setdiff(dfp %>% select(xname,interaction_name_2),dfn %>% select(xname,interaction_name_2)) %>% pull(interaction_name_2)

df1 <- plotDataPrepare(PBMC_normal,1:2,c(1:6,9,10,11:16)) %>% mutate(conditions="normal")
df4 <- plotDataPrepare(PBMC_severe,7:8,c(1:6,9,10,11:16)) %>% mutate(conditions="severe")

df1_source <- plotDataPrepare(PBMC_normal,2,c(3:8,10:15))
df1_target <- plotDataPrepare(PBMC_normal,c(3:8,10:15),2)
df1 <- bind_rows(df1_source,df1_target)
df1$source.target2  <- paste(df1$source, df1$target, sep = " - ") 
df1 <- df1 %>% mutate(conditions="normal")
df1 <- df1 %>% mutate(xname = ifelse(source == "CCR1+ CD16+ monocyte",as.character(target),as.character(source)))
df1 <- df1 %>% distinct(xname,interaction_name_2,.keep_all = TRUE)

df2_source <- plotDataPrepare(PBMC_severe,2,c(3:8,10:15))
df2_target <- plotDataPrepare(PBMC_severe,c(3:8,10:15),2)
df2 <- bind_rows(df2_source,df2_target)
df2$source.target2 <- paste(df2$source, df2$target, sep = " - ") 
df2 <- df2 %>% mutate(conditions="severe")
df2 <- df2 %>% mutate(xname = ifelse(source == "CCR1+ CD16+ monocyte",as.character(target),as.character(source)))
df2 <- df2 %>% distinct(xname,interaction_name_2,.keep_all = TRUE)
df<- bind_rows(df1,df2)
dfp <- df %>% filter(conditions=="severe")
dfn <- df %>% filter(conditions=="normal")
select_row1 <- dplyr::setdiff(dfp %>% select(xname,interaction_name_2),dfn %>% select(xname,interaction_name_2)) %>% pull(interaction_name_2)
