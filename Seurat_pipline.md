# Seurat pipline

## 1. Integration Seurat object

SCT version

```R
set.seed(12345)
library("ggplot2")
library("SeuratObject")
library("dplyr")
#library("harmony")
library("Seurat")
#library("future")
library(tidyverse)
library(Matrix)
#library("clustree")

draw_distribution <- function(seurat_obj, outfname) {
    dat.all=colnames(seurat_obj) %>% as.data.frame() %>% 
			cbind(row=seurat_obj$y) %>%
			cbind(col=seurat_obj$x) %>%
			cbind(slice=seurat_obj$tissue)
    p <- ggplot(dat.all,aes(x=as.numeric(dat.all$row),y=as.numeric(dat.all$col))) + geom_point(size=0.5,shape=16) + 
    theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = '#ffffff')) + 
    coord_fixed(ratio=1)

    ggsave(file=outfname, plot=p, width=50, height=50, limitsize = FALSE)
}

args <- commandArgs(T)

in_dir <- args[1] # gem matrix
n_fnumber <- as.integer(args[2]) # sample id
outp <- args[3] # output path
out_plot <- args[4] # QC output path

# coordinate data
inp = read.table('tissue_info.txt',header=F)
names(inp)=c("slice","stage","y","x","reverse")
obj.list2=list()

for(i in 1:4){
#for(i in 1:3){
    sample=inp[i,]
    dat <- readRDS(paste0(
     in_dir,inp$slice[i],"_spatialObj.rds")
    )

    dat <- RenameCells(object = dat, add.cell.id = dat$tissue)
    if(inp$reverse[i]==0){
        dat$y=dat@images[[inp$slice[i]]]@coordinates$row
        dat$x=dat@images[[inp$slice[i]]]@coordinates$col
    }else{
        temp = dat@images[[inp$slice[i]]]@coordinates$row
        dat@images[[inp$slice[i]]]@coordinates$row = dat@images[[inp$slice[i]]]@coordinates$col
        dat@images[[inp$slice[i]]]@coordinates$col = temp
        dat$y=dat@images[[inp$slice[i]]]@coordinates$row
        dat$x=dat@images[[inp$slice[i]]]@coordinates$col
    }
    obj.list2[[i]] <- SCTransform(dat, assay = "Spatial", vst.flavor = "v2")
}

features <- SelectIntegrationFeatures(object.list = obj.list2, nfeatures = n_fnumber)
print(length(features))
sample.list <- PrepSCTIntegration(object.list = obj.list2, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT", anchor.features = features)
combined<- IntegrateData(anchorset = anchors, normalization.method = "SCT")

saveRDS(combined, file=outp)

if (!dir.exists(out_plot)){
    dir.create(out_plot)
}

# cell distribution after filted
draw_distribution(combined, 
                  paste0(out_plot, "01_cell_distributin_before_filted.pdf"))
# vln plot after filted
p2=VlnPlot(combined, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2)
ggsave(file=paste0(out_plot, "02_Visualize_QC_metrics_before_filted.pdf"), 
       plot=p2)
p3=FeatureScatter(combined, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
ggsave(file=paste0(out_plot, "03_Counts_Features_before_filted.pdf"), plot=p3)
```

Normalized

```R
set.seed(12345)
library("ggplot2")
library("SeuratObject")
library("dplyr")
#library("harmony")
library("Seurat")
#library("future")
library(tidyverse)
library(Matrix)
#library("clustree")

draw_distribution <- function(seurat_obj, outfname) {
    dat.all=colnames(seurat_obj) %>% as.data.frame() %>% 
			cbind(row=seurat_obj$y) %>%
			cbind(col=seurat_obj$x) %>%
			cbind(slice=seurat_obj$tissue)
    p <- ggplot(dat.all,aes(x=as.numeric(dat.all$row),y=as.numeric(dat.all$col))) + geom_point(size=0.5,shape=16) + 
    theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = '#ffffff')) + 
    coord_fixed(ratio=1)

    ggsave(file=outfname, plot=p, width=50, height=50, limitsize = FALSE)
}

args <- commandArgs(T)

in_dir <- args[1] # gem matrix
n_fnumber <- as.integer(args[2]) # sample id
outp <- args[3] # output path
out_plot <- args[4] # QC output path

# coordinate data
inp = read.table('tissue_info.txt',header=F)
names(inp)=c("slice","stage","y","x","reverse")
obj.list2=list()

for(i in 1:4){
#for(i in 1:3){
    sample=inp[i,]
    dat <- readRDS(paste0(
     in_dir,inp$slice[i],"_spatialObj.rds")
    )

    dat <- RenameCells(object = dat, add.cell.id = dat$tissue)
    if(inp$reverse[i]==0){
        dat$y=dat@images[[inp$slice[i]]]@coordinates$row
        dat$x=dat@images[[inp$slice[i]]]@coordinates$col
    }else{
        temp = dat@images[[inp$slice[i]]]@coordinates$row
        dat@images[[inp$slice[i]]]@coordinates$row = dat@images[[inp$slice[i]]]@coordinates$col
        dat@images[[inp$slice[i]]]@coordinates$col = temp
        dat$y=dat@images[[inp$slice[i]]]@coordinates$row
        dat$x=dat@images[[inp$slice[i]]]@coordinates$col
    }
}

var.features <- SelectIntegrationFeatures(object.list = obj.list2, nfeatures = 3000)
combined <- reduce(obj.list2,merge)

VariableFeatures(combined) <- var.features
combined <- NormlizeData(combined)
combined <- ScaleData(combined)
saveRDS(combined, file=outp)

if !dir.exists(out_plot){
    dir.create(out_plot)
}

# cell distribution after filted
draw_distribution(combined, 
                  paste0(out_plot, "01_cell_distributin_before_filted.pdf"))
# vln plot after filted
p2=VlnPlot(combined, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2)
ggsave(file=paste0(out_plot, "02_Visualize_QC_metrics_before_filted.pdf"), 
       plot=p2)
p3=FeatureScatter(combined, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
ggsave(file=paste0(out_plot, "03_Counts_Features_before_filted.pdf"), plot=p3)
```

## 2. Prepare cluster (QC & PCA & UMAP)

Different QC level and downstream pipeline

```R
set.seed(12345)
library("ggplot2")
library("SeuratObject")
library("dplyr")
library("harmony")
library("Seurat")
#library("future")
library(tidyverse)
library(Matrix)
#library("clustree")

draw_distribution <- function(seurat_obj, outfname) {
    dat.all=colnames(seurat_obj) %>% as.data.frame() %>% 
			cbind(row=seurat_obj$y) %>%
			cbind(col=seurat_obj$x) %>%
			cbind(slice=seurat_obj$tissue)
    p <- ggplot(dat.all,aes(x=as.numeric(dat.all$row),y=as.numeric(dat.all$col))) + geom_point(size=0.5,shape=16) + 
    theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = '#ffffff')) + 
    coord_fixed(ratio=1)

    ggsave(file=outfname, plot=p, width=50, height=50, limitsize = FALSE)
}

args <- commandArgs(T)

in_obj <- args[1] # input object
nf_min <- as.integer(args[2]) # min nfeature
nf_max <- as.integer(args[3]) # max nfeature
out_path <- paste0(args[4], "/") # output path

combined <- readRDS(in_obj)

if (!dir.exists(out_path)){
    dir.create(out_path)
}

# filting
combined <- subset(combined, subset = nFeature_Spatial > nf_min & nFeature_Spatial < nf_max)

# cell distribution after filted
draw_distribution(combined, 
                  paste0(out_path, "01_cell_distributin_after_filted.pdf"))
# vln plot after filted
p2=VlnPlot(combined, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2)
ggsave(file=paste0(out_path, "02_Visualize_QC_metrics_after_filted.pdf"), 
       plot=p2)
p3=FeatureScatter(combined, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
ggsave(file=paste0(out_path, "03_Counts_Features_after_filted.pdf"), plot=p3)

combined <- RunPCA(combined, verbose = FALSE)
p4 <- ElbowPlot(combined, ndims = 50)
ggsave(file=paste0(out_path, "04_ElbowPlot.pdf"), plot=p4)
combined <- RunHarmony(combined, assay.use="SCT", group.by.vars = "tissue") # 去除批次效应
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
saveRDS(combined, file=paste0(out_path, "prep_clt.rds"))
```

## 3. Clustering

```R
set.seed(12345)
library("ggplot2")
library("SeuratObject")
library("dplyr")
#library("harmony")
library("Seurat")
library("future")
library(tidyverse)
library(Matrix)
library("clustree")

my_color <- c("#F7E1A0","#C075A6","#1c04ba","#AAF400","#1CFFCE",
              "#3283FE","#F8A19F","#66B0FF","#7ED7D1","#90AD1C",
              "#2ED9FF","#DEA0FD","#AA0DFE","#F6222E","#325A9B",
              "#C4451C","#1C8356","#85660D","#B10DA1","#FBE426",
              "#1CBE4F","#FA0087","#FC1CBF","#5A5156","#BDCDFF",
              "#782AB6","#B5EFB5","#FE00FA","#822E1C","#FEAF16",
              "#FF8247","#1C7F93","#D85FF7","#683B79","#B00068",
              "#FF3030","#d14dac","#2dfab9","#FFFF00","#1CE6FF",
              "#FF34FF","#FF4A46","#008941","#006FA6","#A30059",
              "#FFE4E1","#0000A6","#B79762","#004D43","#8FB0FF",
              "#997D87","#5A0007","#809693","#1B4400","#4FC601",
              "#3B5DFF","#BA0900","#FF2F80","#6B7900","#00C2A0",
              "#FFAA92","#FF90C9","#B903AA","#DDEFFF","#7B4F4B",
              "#A1C299","#0AA6D8","#00A087FF","#4DBBD5FF","#E64B35FF",
              "#3C5488FF","#0067A5","#63FFAC","#F38400","#A1CAF1",
              "#C2B280","#848482","#E68FAC", "#F99379", "#604E97",
              "#F6A600", "#B3446C","#DCD300","#882D17", "#8DB600","#654522", 
              "#E25822", "#2B3D26", "#191970","#000080","#6495ED","#1E90FF",
              "#00BFFF","#00FFFF","#FF1493", "#FF00FF","#A020F0","#63B8FF",
              "#008B8B","#54FF9F", "#00FF00","#76EE00","#FFF68F","Yellow1","Gold1",
              "DarkGoldenrod4","#FF6A6A","#FF8247","#FFA54F","#FF7F24", "#FF3030")

plot_umap <- function(df, outprefix, col_mapping){
    p <- ggplot() + 
    geom_point(df, mapping=aes(x=UMAP_1, y=UMAP_2, color=cluster), shape=16, size=1.5) +
    theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = '#ffffff')) +
    scale_colour_manual(values = col_mapping)

    ggsave(outprefix, plot=p, width = 10, height = 10)
}

plot_spatial <- function(df, outprefix, col_mapping){
    p <- ggplot() +
    geom_point(data=df,mapping=aes(x=x,y=y,color=cluster), size=0.5, shape=16) + 
    #scale_fill_gradientn(colours =c("#f6fd2d", "#ff0000")) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_rect(fill="white"),plot.background=element_rect(fill="white")) +
    scale_color_manual(values = col_mapping) + 
    #scale_color_gradientn(colours =c("#555B73","#00F2FF","#1CFF00","#FFF500","#FF4700")) +
    #scale_color_gradientn(fill = Greens) +
    #scale_color_gradientn(colours =c("#fafcc8", "#ff0000")) +
    #scale_fill_distiller(palette = "Greens") +
    #scale_color_gradientn(colours =c("#241333","#312547","#8B1C1D","#FAEA00")) +
    coord_fixed(ratio=1)

  ggsave(outprefix, p, width=50, heigh=50, limitsize = FALSE)
}

args <- commandArgs(T)

in_obj <- args[1] # input object
out_path <- paste0(args[2], "/") # output path

if (!dir.exists(out_path)){
    dir.create(out_path)
}
combined <- readRDS(in_obj)

combined <- FindClusters(combined, resolution = seq(0.2,2,by=0.2), verbose = FALSE)
saveRDS(combined, paste0(out_path, "/clustered.rds"))
#head(combined@meta.data)
# selected columns with res
metadat <- combined@meta.data
umap <- as.data.frame(combined@reductions$umap@cell.embeddings)
umap <- umap[rownames(metadat),]
metadat$UMAP_1 <- umap$UMAP_1
metadat$UMAP_2 <- umap$UMAP_2
meta <- metadat %>% select(contains("_res."))
prefix_clt <- colnames(meta) %>% grep("0.2", ., value = T) %>% gsub("0.2", "", .)
p_tree <- clustree(meta, prefix = prefix_clt)
ggsave(filename = paste0(out_path, "clustree.pdf"), p_tree, width = 20, height = 30)

for (res in c(seq(0.2, 2, 0.2))){

    df <- metadat[,c("UMAP_1", "UMAP_2", "x", "y", paste0(prefix_clt, as.character(res)))]
    colnames(df) <- c("UMAP_1", "UMAP_2", "x", "y", "cluster")
    df$cluster <- as.character(df$cluster)
    colors <- my_color[1:length(names(table(df$cluster)))]
    col_mapping <- setNames(colors, names(table(df$cluster)))

    plot_df <- df[, c("UMAP_1", "UMAP_2", "cluster")]
    plot_umap(plot_df, paste0(out_path, "UMAP.res.", as.character(res), ".pdf"), col_mapping)
    
    plot_df2 <- df[,c("x", "y", "cluster")]
    plot_spatial(plot_df2, paste0(out_path, "spatial.res.", as.character(res), ".pdf"),
                 col_mapping)
}

df3 <- metadat[,c("UMAP_1", "UMAP_2", "tissue")]
colnames(df3) <- c("UMAP_1", "UMAP_2", "cluster")
tissues <- names(table(df3$cluster))
colors <- my_color[1:length(tissues)]
col_mapping <- setNames(colors, tissues)
plot_umap(df3, paste0(out_path, "UMAP.tissue.pdf"), col_mapping)
```

