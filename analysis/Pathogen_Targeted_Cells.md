# Pathogen Targeted Cells

## Distribution of pathogen reads

```R
set.seed(12345)
library("ggplot2")
library("dplyr")
library("Seurat")
library(tidyverse)
library(viridis)

pi_merged <- readRDS("pi_merged.rds")
pi_info <- read.table("inf_pi.tab", header = T, sep="\t")
pi_info_sub <- subset(pi_info, subset=nfeature_10.12=="Y")
pi_merged$cellid <- rownames(pi_merged@meta.data)
pi.metadat <- pi_merged@meta.data
pi.metadat$nfeature10.12 <- ifelse(pi.metadat$cellid %in% pi_info_sub$cellid, "Y", "N")
pi.metadat.sub <- subset(pi.metadat, subset=nfeature10.12=="Y")[,c("cellid", "nfeature10.12", "nCount_Spatial")]
colnames(pi.metadat.sub) <- c("cellid", "pi.nfeature10_12", "pi_Count_Spatial")

potato.merged <- readRDS("final_used_data.RDS")
potato.metadat <- potato.merged@meta.data

merged.dat <- merge(potato.metadat, pi.metadat.sub, by = "cellid", all.x = TRUE)
merged.dat$pi.nfeature10_12 <- ifelse(is.na(merged.dat$pi.nfeature10_12 ), "N", merged.dat$pi.nfeature10_12)
merged.dat$pi_Count_Spatial <- ifelse(is.na(merged.dat$pi_Count_Spatial), 0, merged.dat$pi_Count_Spatial)

dat.plot = merged.dat

cutoff = c(0, 25, 50, 75, 100, 250, 500)
dat.plot$cutoff <- ifelse(dat.plot$pi_Count_Spatial == 0, "l1", "N")
dat.plot$cutoff <- ifelse(dat.plot$pi_Count_Spatial > 500, "l8", dat.plot$cutoff)
for (i in 1:6){
    dat.plot$cutoff <- ifelse((dat.plot$pi_Count_Spatial > cutoff[i] & dat.plot$pi_Count_Spatial <= cutoff[i+1]),
    paste0("l", as.character(i+1)), dat.plot$cutoff)
}

p <- ggplot() +
    geom_point(data=dat.plot, mapping=aes(x=x,y=y, color=cutoff), size=0.5, shape=16) +
    # geom_point(data=dat.plot[dat.plot$pi_Feature_Spatial==0,],mapping=aes(x=row,y=col), color="#b4e9ff", size=0.5, shape=16) + 
    # scale_color_viridis_c(option = "B",direction = 1) +
    scale_color_manual(values =c(viridis(8))) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_rect(fill="black"),plot.background=element_rect(fill="black")) +
    #scale_color_gradientn(colours =c("#555B73","#00F2FF","#1CFF00","#FFF500","#FF4700")) +
    #scale_color_gradientn(fill = Greens) +
    #scale_color_gradientn(colours =c("#fafcc8", "#ff0000")) +
    #scale_fill_distiller(palette = "Greens") +
    #scale_color_gradientn(colours =c("#241333","#312547","#8B1C1D","#FAEA00")) +
    coord_fixed(ratio=1)

ggsave("fig2A.pi_distribution.pdf", width=50, height=50, limitsize=F)

# reads number
# l1 0
# l2 0~25
# l3 25~50
# l4 50~75
# l5 75~100
# l6 100~250
# l7 250~500
# l8 >500
```

Pathogen reads change over time point

```R
# add ck
library(Seurat)
library(ggplot2)
library(tidyverse)

library(patchwork)

time1 <- c("3h"="#ffdf91", "6h"="#eaac7f", "12h"="#91684a", "24h"="#493323")
time1_ck <- c("3h"="#e5f1e3", "6h"="#a3cd9e", "12h"="#529471", "24h"="#35635b")
time2 <-  c("3T"="#ffdf91", "6h"="#eaac7f", "12T"="#91684a", "24T"="#493323")
time2_ck <- c("3T"="#e5f1e3", "6h"="#a3cd9e", "12T"="#529471", "24T"="#35635b")

pi_merged <- readRDS("pi_merged.rds")
pi_info <- read.table("inf_pi.tab", header=T, sep="\t")
pi_info_sub <- subset(pi_info, subset=nfeature_10.12=="Y")
pi_merged$cellid <- rownames(pi_merged@meta.data)
pi.metadat <- pi_merged@meta.data
pi.metadat$nfeature10.12 <- ifelse(pi.metadat$cellid %in% pi_info_sub$cellid, "Y", "N")
pi.metadat.sub <- subset(pi.metadat, subset=nfeature10.12=="Y")[,c("cellid", "nfeature10.12", "nCount_Spatial", "time")]
colnames(pi.metadat.sub) <- c("cellid", "pi.nfeature10_12", "pi_Count_Spatial", "time")
pi.metadat.sub$time <- factor(pi.metadat.sub$time, levels = c("3h","6h","12h","24h"))

potato_merged=readRDS('final_used_data.RDS')
potato.metadat <- potato_merged@meta.data
merged.dat <- left_join(potato.metadat, pi.metadat.sub, by = "cellid")
merged.dat$pi.nfeature10_12 <- ifelse(is.na(merged.dat$pi.nfeature10_12), "N", merged.dat$pi.nfeature10_12)
merged.dat$pi_Count_Spatial <- ifelse(is.na(merged.dat$pi_Count_Spatial), 0, merged.dat$pi_Count_Spatial)
merged.dat$time <- merged.dat$hours

merged.dat.inf <- merged.dat[merged.dat$inf_stat=="Inf", ]
time_tab <- as.data.frame(table(merged.dat.inf$time, merged.dat.inf$pi.nfeature10_12))
colnames(time_tab) <-  c("time", "PI", "cell_n")
time_tab$prop <- time_tab$cell_n / tapply(time_tab$cell_n, time_tab$time, sum)[time_tab$time]
time_tab$time <- factor(time_tab$time, levels = c("3T","6h","12T","24T"))
time_tab_ck <- data.frame(
  time=factor(c("3T","6h","12T","24T"), levels = c("3T","6h","12T","24T")),
  prop=rep(0, times=4)
)

p1 <- ggplot()
p1 <- p1 + geom_line(data = time_tab[time_tab$PI=="Y",], mapping = aes(x=time, y=prop, group=1), size = 1.5, color="black") +
      geom_point(data = time_tab[time_tab$PI=="Y",], mapping = aes(x=time, y=prop, color = time), size=5, shape=18) +
      scale_color_manual(values = time2)
p1 <- p1 + geom_line(data = time_tab_ck, mapping = aes(x=time, y=prop, group=1), size = 1.5, color="black") +
      geom_point(data = time_tab_ck, mapping = aes(x=time, y=prop, color = time), size=5, shape=18) +
      scale_color_manual(values = time2_ck)
p1 <- p1 + theme_bw()

p2 <- ggplot() + 
geom_jitter(data=pi.metadat.sub, mapping = aes(x=time, y=pi_Count_Spatial, color=time), shape=16, size=1.5) + 
scale_color_manual(values = time1)+
theme_bw()

p <- p1 + p2
ggsave(filename = "pi_cell.time_distribution.pdf", p, width = 10, height = 4)
```

## Trajectory analysis for PTCs

```R
set.seed(12345)
library("ggplot2")
library("dplyr")
library("Seurat")
library(tidyverse)
library(viridis)
library(patchwork)
library(monocle)

pi_info <- read.table("inf_pi.tab", header = T, sep="\t")
pi_neighbor <- read.table("pi.neighbor.giotto.tab", header=T, sep="\t")
pi_cell <- pi_info[pi_info$nfeature_10.12=="Y",][["cellid"]]
pi_neighbor_cell <- pi_neighbor[pi_neighbor$nb_cells!="others",][["cell_ID"]]
pi_neighbor_cell <- setdiff(pi_neighbor_cell, pi_cell)

potato.merged <- readRDS("final_used_data.RDS")
potato.merged$pi_nb <- ifelse(potato.merged$cellid %in% pi_cell, "pi", "N")
potato.merged$pi_nb <- ifelse(potato.merged$cellid %in% pi_neighbor_cell, "neighbor", potato.merged$pi_nb)
potato.merged$pi_nb <- ifelse(potato.merged$inf_stat=="CK", "CK", potato.merged$pi_nb)

dat.part <- subset(potato.merged, subset=(pi_nb=="pi")|(pi_nb=="CK" & hours=="3T"))
dat.part$hours2 <- ifelse(dat.part$pi_nb=="CK", "0T", dat.part$hours)
table(dat.part$pi_nb, dat.part$hours2)
Idents(dat.part) <- dat.part$hours2
run_monocle2=function(sbj){
    expr_matrix = as(as.matrix(sbj@assays$Spatial@counts),'sparseMatrix')
    p_data=sbj@meta.data
    p_data$celltype=sbj@active.ident
    f_data=data.frame(gene_short_name=row.names(sbj),row.names=row.names(sbj))
    pd=new('AnnotatedDataFrame',data=p_data)
    fd=new('AnnotatedDataFrame',data=f_data)
    cds=newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.5,
                    expressionFamily=negbinomial.size())

    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    return(cds)
}

cds <- run_monocle2(dat.part)

GM_state <- function(cds,root){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State,pData(cds)$hours2)[,root]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}

# specific input gene list via genelist
run_monocle=function(sbj,genelist, cds,out_name,root){
    sbj=FindVariableFeatures(sbj)
    var.genes <- VariableFeatures(sbj)
    cds <- setOrderingFilter(cds, genelist)
    cds=reduceDimension(cds,max_components=2,method="DDRTree")
    cds=orderCells(cds)
    root_in=GM_state(cds, root)
    cds=orderCells(cds,root_state=root_in)

    p=plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE)
    ggsave(filename=paste0(out_name,".monocle.pseudotime.seuratHVG.pdf"),width=14,height=14)

    p=plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=TRUE)
    ggsave(filename=paste0(out_name,".monocle.state.seuratHVG.pdf"),width=14,height=14)

    p=plot_cell_trajectory(cds,color_by="celltype",size=1,show_backbone=TRUE)
    ggsave(filename=paste0(out_name,".monocle.celltype.seuratHVG.pdf"),width=14,height=14)
    p=plot_cell_trajectory(cds,color_by="hours",size=1,show_backbone=TRUE)
    ggsave(filename=paste0(out_name,".monocle.hours.seuratHVG.pdf"),width=14,height=14)
    
    cds
}

outprefix <- "./"
cds <- run_monocle(dat.part, cds, genelist, paste0(outprefix, 'pi_time'), '0T')

saveRDS(cds, "pi_time.rds")
```

## DGE analysis

DEG identification

```R
# DEG analysis
library(Seurat)
library(ggplot2)
set.seed(12345)
library(viridis)
library(tidyverse)

pi_info <- read.table("inf_pi.tab", header = T, sep="\t")
pi_neighbor <- read.table("pi.neighbor.giotto.tab", header=T, sep="\t")
pi_cell <- pi_info[pi_info$nfeature_10.12=="Y",][["cellid"]]
pi_neighbor_cell <- pi_neighbor[pi_neighbor$nb_cells!="others",][["cell_ID"]]
pi_neighbor_cell <- setdiff(pi_neighbor_cell, pi_cell)

potato.merged <- readRDS("/public/home/daijichen/01_workdir/02_liyuying_2/final_used_data_0807.merged.RDS")
potato.merged$pi_nb <- ifelse(potato.merged$cellid %in% pi_cell, "pi", "N")
potato.merged$pi_nb <- ifelse(potato.merged$cellid %in% pi_neighbor_cell, "pi_nb", potato.merged$pi_nb)
potato.merged$pi_nb <- ifelse(potato.merged$inf_stat=="Noninf", "Noninf", potato.merged$pi_nb)
potato.merged$pi_nb <- ifelse(potato.merged$inf_stat=="CK", "CK", potato.merged$pi_nb)

# pi vs other
potato.merged$pi_nb2 <- ifelse(potato.merged$pi_nb=="pi_nb" | potato.merged$pi_nb=="N" , "other", potato.merged$pi_nb)
table(potato.merged$pi_nb)
table(potato.merged$pi_nb2)

Idents(potato.merged) <- potato.merged$pi_nb2
find_mark <- function(obj, outprefix){
    diff.wilcox = FindAllMarkers(obj, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
    all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
    all.markers$gene=paste0('A157_',all.markers$gene)
    write.table(all.markers,
    file=paste0(outprefix, ".marker_gene.txt"),
    sep="\t",col.names=T,quote=F,row.names=F)
}
diff_exp_paired <- function(obj, id1, id2, outprefix){
    diff.wilcox = FindMarkers(obj, only.pos = FALSE, 
                                min.pct = 0.1, logfc.threshold = 0.1,
                                ident.1 = id1, ident.2 = id2
                                )
    diff.wilcox$gene = rownames(diff.wilcox)
    all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
    all.markers$gene=paste0('A157_',all.markers$gene)
    write.table(all.markers,
    file=paste0(outprefix, ".DEG_gene.txt"),
    sep="\t",col.names=T,quote=F,row.names=F)
}

map(
    c("3T", "6h", "12T", "24T"), 
    function(x){
        sub_obj <- subset(potato.merged, subset=hours==as.character(x))
        find_mark(sub_obj, paste0("inf_pi.pi_other.",x))
        diff_exp_paired(sub_obj, "other", "CK", paste0("inf_pi.pi_other.", x, "other_vs_CK"))
        diff_exp_paired(sub_obj, "pi", "CK", paste0("inf_pi.pi_other.", x, "pi_vs_CK"))
    }
)
```

Pseudotime heatmap plot

```R
# pseudo-time heatmap plot
set.seed(12345)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(monocle)
library(Seurat)
library(viridis)
library(RColorBrewer)

cds <- readRDS("pi_time.seuratHVG.rds")

gene_list <- readRDS("inf_pi.pi_other.time.vs_CK.DEG_gene.rds")
up_s_gene_list <- c()
for (n in grep("up_s", names(gene_list), value = TRUE)){
  tmp_df <- gene_list[[n]]
  tmp_df <- tmp_df[tmp_df$pi==1,]
  up_s_gene_list <- c(up_s_gene_list, tmp_df[["gene"]])
}

DEG.genes <- c(up_s_gene_list) %>% unique() %>% gsub("A157_", "",.)

p <- plot_pseudotime_heatmap(cds[DEG.genes,], 
num_clusters = 4, 
cores = 10, 
show_rownames = F, 
return_heatmap=T,
hmcols = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(50))
)

p_order <- p$tree_row$order
p_label <- p$tree_row$label
p_label <- p_label[p_order] # reorder the label
genecluster <- as.data.frame(cutree(p$tree_row, k = 4))
colnames(genecluster) <- c("cluster")
genecluster$gene <- rownames(genecluster)
genecluster <- genecluster[p_label,]
write.table(genecluster, file="gene_clusters.pesudo_order_20240427.txt",
col.names = T, row.names = F, sep = "\t", quote = F)
pdf("inf_pi.heatmap.pesudo_order.20240427.pdf", width = 5, height = 7.5)
print(p)
dev.off()

newdata <- seq(min(pData(cds)$Pseudotime), max(pData(cds)$Pseudotime),length.out = 100)
cell_pseudo_df <- data.frame(Pseudotime = pData(cds)$Pseudotime, celltype = pData(cds)$celltype)

sep <- (newdata[2] - newdata[1])/2
plot_df = cell_pseudo_df[(cell_pseudo_df$Pseudotime >= newdata[1] & cell_pseudo_df$Pseudotime < newdata[1] + sep ), ]
plot_df = as.data.frame(table(plot_df$celltype))
# plot_df$celltype = factor(plot_df$celltype) 必须为因子类型
plot_df$group <- "box1"

for (i in seq(2,100)){
    if (i == 100){
        tmp_df = cell_pseudo_df[((cell_pseudo_df$Pseudotime >= newdata[i]-sep) & (cell_pseudo_df$Pseudotime <= newdata[i])), ]
        tmp_df = as.data.frame(table(tmp_df$celltype))
        tmp_df$group <- paste0("box", as.character(i))
        plot_df <- rbind(plot_df, tmp_df)
    }else {
        tmp_df = cell_pseudo_df[((cell_pseudo_df$Pseudotime >= newdata[i]-sep) & (cell_pseudo_df$Pseudotime < newdata[i] + sep)), ]
        tmp_df = as.data.frame(table(tmp_df$celltype))
        tmp_df$group <- paste0("box", as.character(i))
        plot_df <- rbind(plot_df, tmp_df)
    }
}

# 山脊图
library(ggridges)
library(tidyverse)

df_list <- list()
for (clt in names(table(pData(cds)$celltype))){
    tmp_df <- plot_df[plot_df$Var1==clt, ]
    tmp_df <- tmp_df %>% mutate(expanded = map2(group, Freq, rep)) %>%
    unnest(expanded)
    tmp_df <- as.data.frame(tmp_df)
    tmp_df$expanded <- as.integer(gsub("box", "", tmp_df$expanded))
    df_list[[clt]] <- tmp_df
}
combined_df <- do.call(rbind, df_list)

p <- ggplot(combined_df, aes(x=expanded, y = Var1, fill = Var1)) +
geom_density_ridges() +
theme_ridges() +
scale_fill_manual(values = c("0T"="#1CBE4F", "3T"="#1e90ff", "6h"="#008b8b", "12T"="#326a9b", "24T"="#364958")) +
theme(legend.position = "none") + scale_x_continuous(breaks=seq(0, 100, 20))

ggsave(file="ggridges.pesudo_order.20240427.pdf", plot=p, width = 5, height = 1.5)
```

## Local Spatial Gradient Inference (LSGI) analysis

```R
# we recommend using NMF as a default summarization method
# for gene expression matrix in this example we use RcppML
# for NMF, other methods are acceptable
library(Seurat)
library(Matrix)
library(RcppML)  # https://github.com/zdebruine/RcppML
library(ggplot2)
library(dplyr)
library(LSGI)
library(viridis)
library(RColorBrewer)
set.seed(12345)

# prepare data
potato_merged <- readRDS("final_used_data.RDS")

data <- subset(potato_merged, subset=(hours=="12T") & (inf_stat=="Inf"))
# data <- subset(potato_merged, subset=(inf_stat=="Inf"))
# dim(data@meta.data)

# random remove cells 讨厌质数
celldf <- data.frame(cell=rownames(data@meta.data))
remove_cells <- sample_n(celldf, 67)[["cell"]] # 8167 - 67 = 8100 # all 16
data <- subset(data, subset=(cellid %in% remove_cells), invert=T)

# how much nmf
scan.nmf.mse <- function(obj, ranks = seq(1, 20, 2), tol = 1e-04) {
    # users can customize the scan by changing 'ranks'
    dat <- as.matrix(obj@assays$Spatial@data)
    errors <- c()
    ranks <- seq(1, 30, 2)
    for (i in ranks) {
        # cat('rank: ', i, '\n')
        mod <- RcppML::nmf(data=dat, k=i, tol = 1e-04, verbose = F)
        mse_i <- mse(dat, w=mod$w, d=mod$d, h=mod$h)
        errors <- c(errors, mse_i)
    }
    results <- data.frame(rank = ranks, MSE = errors)
    return(results)
}

scan.nmf.res <- scan.nmf.mse(obj = data)
p <- ggplot(scan.nmf.res, aes(x = rank, y = MSE)) + geom_point(size = 0.7) +
    geom_smooth(method = "loess", span = 0.2, color = "black",
        linewidth = 1, se = F) + labs(x = "NMF rank", y = "MSE") +
    theme_classic() + scale_y_continuous(expand = c(0.01, 0)) +
    theme(aspect.ratio = 1)
ggsave("nmf_scan.pdf", width=10, height=10)

# take rank as 9 due to the rank plot
sr.nmf <- function(obj, k = 10, tol = 1e-06, assay = "Spatial") {
    dat <- as.matrix(obj@assays$Spatial@data)
    nmf_model <- RcppML::nmf(data=dat, k=k, tol = 1e-04, verbose = F)
    embeddings <- t(nmf_model$h)
    rownames(embeddings) <- colnames(obj)
    colnames(embeddings) <- paste0("nmf_", 1:k)
    loadings <- nmf_model$w
    rownames(loadings) <- rownames(obj)
    obj@reductions$nmf <- CreateDimReducObject(embeddings = embeddings,
        loadings = loadings, key = "nmf_", assay = assay)
    return(obj)
}

data_nmf <- sr.nmf(obj = data, k = 9, tol = 1e-05)

# input LSGI
spatial_coords <- data.frame(X=data_nmf@meta.data[["y"]], Y= data_nmf@meta.data[["x"]])
embeddings <- data_nmf@reductions$nmf@cell.embeddings
rownames(spatial_coords) <- rownames(embeddings)

lsgi.res <- local.traj.preprocessing(spatial_coords = spatial_coords, embeddings = embeddings)

# visual
p <- plt.factors.gradient.ind(info = lsgi.res, r_squared_thresh = 0.4, minimum.fctr = 10, arrow.length.scale = 0.2)
p <- p + coord_fixed(ratio=1)
ggsave("nmf_landscape.pdf", width=20, height=20, limitsize=F)

# my plot function
plt.factor.gradient.ind_my <- function(info, fctr, r_squared_thresh, arrow.length.scale = 1){

  lin.res.df <- get.ind.rsqrs(info)
  lin.res.df <- na.omit(lin.res.df)
  grid.info <- info$grid.info
  embd <- info[["embeddings"]]
  spatial_coords <- info$spatial_coords
  spatial_coords$factor <- embd[, fctr]
  p <- ggplot(spatial_coords, aes(x=X, y=Y, color = factor)) +
    geom_point(size=2, shape = 20, stroke = 0) +
    scale_color_gradientn(colours = c(colorRampPalette(c("#c9caca", "#fa1e1a"))(30)))+
    geom_segment(data = lin.res.df[lin.res.df$fctr == fctr & lin.res.df$rsquared > r_squared_thresh, ],
                 aes(xend = X + vx.u*arrow.length.scale, yend = Y + vy.u*arrow.length.scale, fill = NULL),
                 color = "black",
                 linewidth = 0.4, arrow = arrow(length = unit(0.1, "cm"))) +
    theme_classic() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_rect(fill="white"),plot.background=element_rect(fill="white")) +
    theme(axis.text.x = element_text(face = "bold", color = "black",
                                     size = 12, angle = 0, hjust = 1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 12, angle = 0))
  return(p)
}

# echo nmf
for (i in 1:9){
    print(i)
    p <- plt.factor.gradient.ind_my(info = lsgi.res, fctr = paste0("nmf_", as.character(i)), r_squared_thresh = 0.4, arrow.length.scale = 0.3)
    p <- p + coord_fixed(ratio=1)
    ggsave(paste0("nmf_", as.character(i), ".pdf"), width=20, height=20, limitsize=F)
}

# feature
get.nmf.info <- function(obj, top.n = 50) {
    feature.loadings <- as.data.frame(obj@reductions$nmf@feature.loadings)

    top.gene.list <- list()
    for (i in 1:ncol(feature.loadings)) {
        o <- order(feature.loadings[, i], decreasing = T)[1:top.n]
        features <- rownames(feature.loadings)[o]
        top.gene.list[[colnames(feature.loadings)[i]]] <- features
    }
    nmf.info <- list(feature.loadings = feature.loadings, top.genes = top.gene.list)
    return(nmf.info)
}

nmf_info <- get.nmf.info(data_nmf, top.n=50)
nmf_list <- list()
for (n in names(nmf_info$top.genes)){
    nmf_list[[n]] <- data.frame(nmf=rep(n,times=50), gene = nmf_info$top.genes[[n]])
}
nmf_top_df = do.call(rbind, nmf_list)
write.table(nmf_top_df, "top_gene.tab", col.names = T, row.names = F, sep = "\t", quote = F)

# nmf matrix
saveRDS(data_nmf@reductions$nmf@cell.embeddings, "nmf_info.rds")
```

