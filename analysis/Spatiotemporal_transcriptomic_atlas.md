# Spatiotemporal transcriptomic atlas

## Base statistic for UMI and gene features

```R
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggdist)

my_color = c(
    "0"="#93d0fc", "1"="#2aaae1", "2"="#3f5246", "3"="#ad3074", "4"="#E7B10A", 
    "5"="#91a3b0", "6"="#d58282", "7"="#418557", "8"="#6dc831", "9"="#eacd01",
    "10"="#ff9f22", "11"="#af5b2b", "12"="#2e5073", "13"="#E6D2AA", "14"="#42032C",
    "15"="#752f9a", "16"="#721c28", "17"="#2c2cff", "18"="#d1bac8", "19"="#009580",
    "20"="#94ac78", "21"="#9f59f8", "22"="#ffff00", "23"="#fc0cd4", "24"="#ffccff"
)
my_color3 = c(
    "Inf" = "#C4451C",
    "Noninf" = "#3283FE",
    "CK" = "#1C8356"
)
infect_region = c("3", "4", "6", "10", "11", "13", "14", "16", "23")
noninfect_region = c("0", "1", "2", "5", "7", "8", "9", "12", "15", "17", "18", "19", "20",
                    "21", "22", "24")
clt_order = c("3", "4", "6", "10", "11", "13", "14", "16", "23",
         "0", "1", "2", "5", "7", "8", "9", "12", "15", "17", "18", "19", "20",
         "21", "22", "24")
time_color <- c("0T"="#1CBE4F", "3T"="#1e90ff", "6h"="#008b8b", "12T"="#326a9b", "24T"="#364958")

potato.merged <- readRDS('final_used_data.RDS')

plot_df <- data.frame(time=potato.merged$hours, nFeature_Spatial=potato.merged$nFeature_Spatial, nCount_Spatial=potato.merged$nCount_Spatial)
plot_df$time <- factor(plot_df$time, levels=c("3T", "6h", "12T", "24T"))

p1 <- ggplot(plot_df, aes(time, nFeature_Spatial)) + 
  ggdist::stat_halfeye(aes(fill=time),adjust = .5, width = .4, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(aes(color=time),width = .2, outlier.shape = 21) + 
#  geom_jitter(aes(color=time),width = .05, alpha = .3) +
  coord_flip()+
  scale_fill_manual(values=time_color) + 
  scale_color_manual(values=time_color) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 
p2 <- ggplot(plot_df, aes(time, nCount_Spatial)) + 
  ggdist::stat_halfeye(aes(fill=time),adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(aes(color=time),width = .2, outlier.shape = 21) + 
#  geom_jitter(aes(color=time),width = .05, alpha = .3) +
  coord_flip()+
  scale_fill_manual(values=time_color) + 
  scale_color_manual(values=time_color) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 

p <- p1 + p2

ggsave(filename="figS_time_boxplot.umi_nfeature.pdf", p, width=12, height=7)
```

## Clusters in UMAP and spatial plot

```R
# UMAP
library(Seurat)
library(ggplot2)

my_color2 = c(
    "#93d0fc", "#2aaae1", "#3f5246", "#ad3074", "#ff4c4c", 
    "#91a3b0", "#d58282", "#418557", "#6dc831", "#eacd01",
    "#ff9f22", "#af5b2b", "#2e5073", "#e06666", "#af2b3d",
    "#752f9a", "#721c28", "#2c2cff", "#d1bac8", "#009580",
    "#94ac78", "#9f59f8", "#ffff00", "#fc0cd4", "#ffccff"
)

my_color = c(
    "0"="#93d0fc", "1"="#2aaae1", "2"="#3f5246", "3"="#ad3074", "4"="#E7B10A", 
    "5"="#91a3b0", "6"="#d58282", "7"="#418557", "8"="#6dc831", "9"="#eacd01",
    "10"="#ff9f22", "11"="#af5b2b", "12"="#2e5073", "13"="#E6D2AA", "14"="#42032C",
    "15"="#752f9a", "16"="#721c28", "17"="#2c2cff", "18"="#d1bac8", "19"="#009580",
    "20"="#94ac78", "21"="#9f59f8", "22"="#ffff00", "23"="#fc0cd4", "24"="#ffccff"
)

potato.merged <- readRDS("final_used_data.RDS")

umap_df <- data.frame(
    cluster = factor(as.character(potato.merged$seurat_clusters), levels = as.character(0:24)),
    UMAP_1 = potato.merged$UMAP_1,
    UMAP_2 = potato.merged$UMAP_2
)

p <- ggplot() +
geom_point(umap_df, mapping=aes(x=UMAP_1, y=UMAP_2, color=cluster), shape=16, size=2) +
scale_color_manual(values = my_color) +
theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = '#ffffff'))

ggsave("fig1D_UMAP.pdf", p, width = 20, height = 15)
```

```R
# high light
library(Seurat)
library(ggplot2)
library(parallel)
plan("multicore", workers = 5)

potato.merged <- readRDS("final_used_data.merged.RDS")

plot_df <- data.frame(
    cluster = factor(as.character(potato.merged$seurat_clusters), levels = as.character(0:24)),
    UMAP_1 = potato.merged$UMAP_1,
    UMAP_2 = potato.merged$UMAP_2,
    hours = potato.merged$hours,
    row=potato.merged$y,
    col=potato.merged$x
)

cluster_name <- names(table(potato.merged$seurat_clusters))
#cluster_name <- c("16", "6")
highlight_spaial <- function(i){
    clt=cluster_name[i]
    p <- ggplot() + 
    geom_point(plot_df[plot_df != clt, ], mapping=aes(x=row, y=col), color="#afd3e0", shape=16, size=0.5) +
    geom_point(plot_df[plot_df == clt, ], mapping=aes(x=row, y=col), color="#ee5a19", shape=16, size=0.5) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    panel.background = element_rect(fill="black"),plot.background=element_rect(fill="black"))

    ggsave(paste0("cluster_highlight.spatial.", "clt", clt, ".pdf"), width=50, height=30, limitsize = FALSE)
}
highlight_umap <- function(i){
    clt=cluster_name[i]
    p <- ggplot() + 
    geom_point(plot_df[plot_df != clt, ], mapping=aes(x=UMAP_1, y=UMAP_2), color="#afd3e0", shape=16, size=2) +
    geom_point(plot_df[plot_df == clt, ], mapping=aes(x=UMAP_1, y=UMAP_2), color="#ee5a19", shape=16, size=2) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    panel.background = element_rect(fill="white"),plot.background=element_rect(fill="white"))

    ggsave(paste0("cluster_highlight.umap", "clt", clt, ".pdf"), width=20, height=20, limitsize = FALSE)
}


mclapply(1:length(cluster_name),highlight_spaial,mc.cores=5)
mclapply(1:length(cluster_name),highlight_umap,mc.cores=5)
```

## Gene expression correlation between clusters

```R
set.seed(12345)
library(ggplot2)
library(Seurat)
library(pheatmap)

library(viridis)

callback = function(hc, mat){
  dend = dendextend::rotate(hc, order = rev(as.character(c(13,16,6,10,3,14,4,11,8,0,1,2,7,5,12,15,9))))
  as.hclust(dend)
}

potato.merged <- readRDS("final_used_data.RDS")
Idents(potato.merged) <- potato.merged$seurat_clusters
clt_remove <- c(seq(17,24))
potato.merged <- subset(potato.merged, subset= seurat_clusters %in% clt_remove, invert=T)
potato.merged <- FindVariableFeatures(potato.merged,nfeatures=2000)
table(potato.merged$seurat_clusters)

var.gene <- VariableFeatures(potato.merged)
obj_exp <- AverageExpression(potato.merged, group.by="seurat_clusters", features = var.gene)
exp.mat <- obj_exp$Spatial

pdf("fig1G.cor.pdf")
pheatmap(cor(exp.mat),treeheight_row = 50,treeheight_col = 15,show_rownames=T,show_colnames=T,
clustering_distance_rows = "correlation",clustering_distance_cols = "correlation", clustering_callback = callback,
color = c("#440154FF", "#440154FF", "#440154FF", "#440154FF", "#440154FF", viridis(30)))

dev.off()

```

## Marker gene visulation

Marker genes were identified with **diff_marker_exp** function in `pipeline/seurat_pipeline.md`

```R
library(Seurat)
library(ggplot2)
library(RColorBrewer)

plot_dot_heatmap <- function(obj, genelist, clusterorder){
    p <- DotPlot(object = obj, features = genelist, idents=clusterorder) + 
    theme_classic() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank()) +
    theme(axis.text.x = element_text(hjust = 1, vjust = .5,angle = 90,size=7),
    panel.background = element_rect(fill="white"),
    axis.title.x = element_text(colour = "white")) +
    #scale_colour_gradient2(low = "#4f79b7", mid = "lightgrey", high = "red") + 
    scale_colour_gradientn(colors=rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(11))) + 
    scale_size(range = c(2,10)) +
    coord_flip()
    ggsave("fig1F_dotheatmap.20240126.pdf", width = 10, height = 15, limitsize = FALSE)
}

potato.merged <- readRDS("final_used_data.RDS")

genelist <- read.table("cluster.marker_gene.txt", header = T, sep = "\t")
genelist$gene <- gsub("A157_", "", genelist$gene)

clt_order <- as.character(c(15,9,13,12, 2, 1, 0, 10,3,4,6,16,5,14,8,7,11))

potato.merged$seurat_clusters <- factor(
as.character(potato.merged$seurat_clusters), levels=c(clt_order, 17:24))
Idents(potato.merged) <- potato.merged$seurat_clusters

plot_dot_heatmap(potato.merged, genelist, clt_order)
```

