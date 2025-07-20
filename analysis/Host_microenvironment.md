# Dynamic responses within host microenvironment

## Trajectory analysis

Monocle2

 ```R
 # monocle2
 library(monocle)
 library(Seurat)
 library(tidyverse)
 library(ggplot2)
 library(viridis)
 args <- commandArgs(T)
 
 my_color = c(
     "0"="#93d0fc", "1"="#2aaae1", "2"="#3f5246", "3"="#ad3074", "4"="#E7B10A", 
     "5"="#91a3b0", "6"="#d58282", "7"="#418557", "8"="#6dc831", "9"="#eacd01",
     "10"="#ff9f22", "11"="#af5b2b", "12"="#2e5073", "13"="#E6D2AA", "14"="#42032C",
     "15"="#752f9a", "16"="#721c28", "17"="#2c2cff", "18"="#d1bac8", "19"="#009580",
     "20"="#94ac78", "21"="#9f59f8", "22"="#ffff00", "23"="#fc0cd4", "24"="#ffccff"
 )
 
 potato_merged=readRDS(args[1])
 outprefix <- args[2]
 
 potato_merged$seurat_clusters <- as.character(potato_merged$seurat_clusters)
 Idents(potato_merged) <- potato_merged$seurat_clusters
 
 if (!dir.exists(outprefix)){
   dir.create(outprefix)
 }
 
 run_monocle2=function(sbj){
     expr_matrix = as(as.matrix(sbj@assays$Spatial@counts),'sparseMatrix')
     p_data=sbj@meta.data
     p_data$celltype=sbj@active.ident # 当前活动的idents changed
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
 
 GM_state <- function(cds,root){
     if (length(unique(pData(cds)$State)) > 1){
         T0_counts <- table(pData(cds)$State,pData(cds)$celltype)[,root]
         return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
     } else {
         return (1)
     }
 }
 
 # specific input gene list
 run_monocle=function(sbj,cds, diff.genes, out_name,root){
     cds <- setOrderingFilter(cds, diff.genes)
     cds=reduceDimension(cds,max_components=2,method="DDRTree")
     cds=orderCells(cds)
     root_in=GM_state(cds, root)
     cds=orderCells(cds,root_state=root_in)
 
     p=plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE,cell_size=2) + scale_color_gradientn(colours =c(viridis(20)))
     ggsave(filename=paste0(out_name,".monocle.pseudotime.seuratDEG.pdf"),width=14,height=14)
     p=plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=TRUE)
     ggsave(filename=paste0(out_name,".monocle.state.seuratDEG.pdf"),width=14,height=14)
     p=plot_cell_trajectory(cds,color_by="celltype",size=1,show_backbone=TRUE) + scale_color_manual(values=my_color)
     ggsave(filename=paste0(out_name,".monocle.celltype.seuratDEG.pdf"),width=14,height=14)
     p=plot_cell_trajectory(cds,color_by="inf_stat",size=1,show_backbone=TRUE)
     ggsave(filename=paste0(out_name,".monocle.inf_stat.seuratDEG.pdf"),width=14,height=14)
     
     cds
 }
 
 # all cells
 print("all cells")
 cds <- run_monocle2(potato_merged)
 
 cds <- run_monocle(potato_merged, cds, paste0(outprefix, "potato_merged.inf2noninf"), "16")
 saveRDS(cds, "potato_merged.monocle2.seuratDEG.pdf.rds")
 ```

Monocle3

```R
# monocle3
library(Seurat)
library(monocle3)
set.seed(12345)

args <- commandArgs(T)
data_set <- args[1]
outprefix <- args[2]
if (!dir.exists(outprefix)){
  dir.create(outprefix)
}

potato_merged=readRDS(data_set)
potato_merged$seurat_clusters <- as.character(potato_merged$seurat_clusters)
Idents(potato_merged) <- potato_merged$seurat_clusters

dat_part <- potato_merged

run_monocle3 <- function(sbj){
    expr_matrix = as(as.matrix(sbj@assays$Spatial@counts),'sparseMatrix')
    p_data=sbj@meta.data
    p_data$celltype=sbj@active.ident # 聚类的结果
    f_data=data.frame(gene_short_name=row.names(sbj),row.names=row.names(sbj))
    pd=new('AnnotatedDataFrame',data=p_data)
    fd=new('AnnotatedDataFrame',data=f_data)
    cds=new_cell_data_set(expr_matrix,
                          cell_metadata = p_data,
                          gene_metadata = f_data)
    return(cds)
}

get_earliest_principal_node <- function(cds, root){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == root) # 使用idents的列
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# all cells
cds <- run_monocle3(dat_part)

reducedDim(cds, type = "PCA") <- dat_part@reductions$pca@cell.embeddings
cds@preprocess_aux$prop_var_expl <- dat_part@reductions$pca@stdev

# umap in cell.embeddings
# cds@int_colData@listData$reducedDims$UMAP <- dat_part@reductions$umap@cell.embeddings
# umap in meta.data
umap.tab <- data.frame(UMAP_1=dat_part$UMAP_1, UMAP_2=dat_part$UMAP_2)
rownames(umap.tab) <- rownames(dat_part@meta.data)
cds@int_colData@listData$reducedDims$UMAP <- umap.tab
cds@clusters$UMAP_so$clusters <- Idents(dat_part)
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)

# must set to NULL so it could continue learn graph 
rownames(cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(cds@int_colData@listData$reducedDims$UMAP) <- NULL
# calculate presudotime
cds <- learn_graph(cds, use_partition = T)
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, '16'))
cds$seurat_clusters <- as.character(cds$seurat_clusters)

write.table(pData(cds), paste0(outprefix, ".monocle3.pData.table"), col.names=T, row.names=F, sep="\t", quote=F)

pdf(paste0(outprefix, ".monocle3.celltype.pdf"))
plot_cells(cds,
           color_cells_by = "inf_stat",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           show_trajectory_graph=F)
dev.off()
pdf(paste0(outprefix, ".monocle3.seurat_clts.pdf"))
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           show_trajectory_graph=F)
dev.off()
pdf(paste0(outprefix, ".monocle3.pseudotime.pdf"))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           show_trajectory_graph=F)
dev.off()
pdf(paste0(outprefix, ".monocle3.state.pdf"))
plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           show_trajectory_graph=F)
dev.off()

saveRDS(cds, file=paste0(outprefix, ".monocle3.rds"))
```

## DEG analysis

```R
# DEG analysis

library(Seurat)
library(tidyverse)
library(ggplot2)

potato_merged=readRDS('final_used_data.RDS')
input_clt <- c(16, 6, 10, 11, 8, 7, 19, 2)
dat_part <- subset(potato_merged, subset=(inf_stat=="CK")|((inf_stat!="CK")&(seurat_clusters %in% input_clt)))
dat_part$DEG <- ifelse(dat_part$inf_stat=="CK", "CK", dat_part$seurat_clusters)
#dat_part$DEG <- ifelse((dat_part$seurat_clusters==16 | dat_part$seurat_clusters==6)&dat_part$inf_stat!="CK", "inf", dat_part$DEG)
#dat_part$DEG <- ifelse((dat_part$seurat_clusters==10 | dat_part$seurat_clusters==11 | dat_part$seurat_clusters==8)&dat_part$inf_stat!="CK", "mid", dat_part$DEG)
#dat_part$DEG <- ifelse((dat_part$seurat_clusters==7 | dat_part$seurat_clusters==19 | dat_part$seurat_clusters==2)&dat_part$inf_stat!="CK", "noninf", dat_part$DEG)
Idents(dat_part) <- dat_part$DEG
#Idents(dat_part) <- dat_part$seurat_clusters

dat.part.3h <- subset(dat_part, subset=(hours=="3T"))
dat.part.6h <- subset(dat_part, subset=(hours=="6h"))
dat.part.12h <- subset(dat_part, subset=(hours=="12T"))
dat.part.24h <- subset(dat_part, subset=(hours=="24T"))

diff_exp_paired <- function(obj, id1, id2, outprefix){
    diff.wilcox = FindMarkers(obj, only.pos = FALSE, 
                                min.pct = 0.25, logfc.threshold = 0.1,
                                ident.1 = id1, ident.2 = id2
                                )
    diff.wilcox$gene = rownames(diff.wilcox)
    all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
    all.markers$gene=paste0('A157_',all.markers$gene)
    write.table(all.markers,
    file=paste0(outprefix, ".DEG_gene.txt"),
    sep="\t",col.names=T,quote=F,row.names=F)
}
diff_exp <- function(obj, outprefix){
    diff.wilcox = FindAllMarkers(obj, only.pos = T, 
                                min.pct = 0.25, logfc.threshold = 0.1
                                )
    all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
    all.markers$gene=paste0('A157_',all.markers$gene)
    write.table(all.markers,
    file=paste0(outprefix, ".marker_gene.txt"),
    sep="\t",col.names=T,quote=F,row.names=F)
}
diff_exp(dat.part.3h, "3h.inf_mid_noninf_cltmrk")
diff_exp(dat.part.6h, "6h.inf_mid_noninf_cltmrk")
diff_exp(dat.part.12h, "12h.inf_mid_noninf_cltmrk")
diff_exp(dat.part.24h, "24h.inf_mid_noninf_cltmrk")
diff_exp(dat_part, "all.inf_mid_noninf_cltmrk")
```

## Pseudotime heatmap plot

```R
# this script used to add cell distribution upon heatmap 
# eg. 12h cluster8 inf noninf cells upon pesudotime split 

library(monocle)
library(Seurat)
library(tidyverse)

my_color = c(
      "0"="#93d0fc", "1"="#2aaae1", "2"="#3f5246", "3"="#ad3074", "4"="#E7B10A", 
      "5"="#91a3b0", "6"="#d58282", "7"="#418557", "8"="#6dc831", "9"="#eacd01",
      "10"="#ff9f22", "11"="#af5b2b", "12"="#2e5073", "13"="#E6D2AA", "14"="#42032C",
      "15"="#752f9a", "16"="#721c28", "17"="#2c2cff", "18"="#d1bac8", "19"="#009580",
      "20"="#94ac78", "21"="#9f59f8", "22"="#ffff00", "23"="#fc0cd4", "24"="#ffccff"
  )

args <- commandArgs(T)

potato_merged=readRDS(args[1])
potato_merged$DEG <- ifelse(potato_merged$seurat_clusters==16 | potato_merged$seurat_clusters==6, "inf", "N")
potato_merged$DEG <- ifelse(potato_merged$seurat_clusters==10 | potato_merged$seurat_clusters==11 | potato_merged$seurat_clusters==8, "mid", potato_merged$DEG)
potato_merged$DEG <- ifelse(potato_merged$seurat_clusters==7 | potato_merged$seurat_clusters==19 | potato_merged$seurat_clusters==2, "noninf", potato_merged$DEG)
Idents(potato_merged) <- potato_merged$DEG
potato_merged$seurat_clusters <- as.character(potato_merged$seurat_clusters)
gene_f <- args[2]

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

cds <- run_monocle2(potato_merged)

GM_state <- function(cds,root){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State,pData(cds)$seurat_clusters)[,root]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}

run_seuratDEG=function(sbj,cds,out_name,root){
    sbj <- NormalizeData(sbj, normalization.method = "LogNormalize", scale.factor = 10000)
    sbj <- ScaleData(sbj) 
    deg.cluster <- FindAllMarkers(sbj)
    diff.genes <- subset(deg.cluster,p_val<0.05)$gene
    cds <- setOrderingFilter(cds, diff.genes)
    cds=reduceDimension(cds,max_components=2,method="DDRTree")
    cds=orderCells(cds)
    root_in=GM_state(cds, root)
    cds=orderCells(cds,root_state=root_in)

    p=plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE,cell_size=2)
    ggsave(filename=paste0(out_name,".monocle.pseudotime.seuratDEG.pdf"),width=14,height=14)
    p=plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=TRUE)
    ggsave(filename=paste0(out_name,".monocle.state.seuratDEG.pdf"),width=14,height=14)
    p=plot_cell_trajectory(cds,color_by="celltype",size=1,show_backbone=TRUE)
    ggsave(filename=paste0(out_name,".monocle.celltype.seuratDEG.pdf"),width=14,height=14)
    p=plot_cell_trajectory(cds,color_by="seurat_clusters",size=1,show_backbone=TRUE)
    ggsave(filename=paste0(out_name,".monocle.seurat_clts.seuratDEG.pdf"),width=14,height=14)

    cds
}

cds <- run_seuratDEG(potato_merged, cds, "inf_mid_noninf", "16")

# immune gene innner time
immune_gene <- read.table(gene_f, sep="\t", header = T)[["gene"]] %>% gsub("A157_", "", .) %>% unique(.)

# plot pesudo time heatmap
p <- plot_pseudotime_heatmap(cds[immune_gene,], 
num_clusters = 10, 
cores = 10, 
show_rownames = F, 
return_heatmap=T)

p_order <- p$tree_row$order
p_label <- p$tree_row$label
p_label <- p_label[p_order] # reorder the label
genecluster <- as.data.frame(cutree(p$tree_row, k = 10))
colnames(genecluster) <- c("cluster")
genecluster$gene <- rownames(genecluster)
genecluster <- genecluster[p_label,]
write.table(genecluster, file=paste0("inf_mid_noninf", ".pesudo_order_20231026.tab"),
col.names = T, row.names = F, sep = "\t", quote = F)

saveRDS(p, file=paste0("inf_mid_noninf", ".pheatmap.pesudo_order_20231026.rds"))
pdf(paste0("inf_mid_noninf", ".pheatmap.pesudo_order_20231026.pdf"), width = 7.5, height = 7.5)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off() 

# prepare data
newdata <- seq(min(pData(cds)$Pseudotime), max(pData(cds)$Pseudotime),length.out = 100)
cell_pseudo_df <- data.frame(Pseudotime = pData(cds)$Pseudotime, celltype = pData(cds)$seurat_clusters)

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

df_list <- list()
for (clt in names(table(plot_df$Var1))){
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
scale_fill_manual(values = my_color) +
theme(legend.position = "none") + scale_x_continuous(breaks=seq(0, 100, 20))

ggsave(file=paste0("inf_mid_noninf", ".ggridges.pesudo_order_20231102.pdf"), plot=p, width = 5, height = 1.5)
```

## Key TF identification

```R

set.seed(12345)
library("dplyr")
library("Seurat")
library("future")
library(tidyverse)
library(Matrix)
options(future.globals.maxSize = 240000 * 1024^2)

# diff exp paired
diff_exp_paired <- function(obj, id1, id2, outprefix){
    diff.wilcox = FindMarkers(obj, only.pos = F, 
                                min.pct = 0.25, logfc.threshold = 0.1,
                                ident.1 = id1, ident.2 = id2
                                )
    diff.wilcox$gene = rownames(diff.wilcox)
    all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
    all.markers$gene=paste0('A157_',all.markers$gene)
    write.table(all.markers,
    file=paste0(outprefix, "_marker_gene.txt"),
    sep="\t",col.names=T,quote=F,row.names=F)
}

args <- commandArgs(T)
clt1 <- args[1]
clt2 <- args[2]

obj <- readRDS("../Dataset2.6h.20231022.rds")
diff_exp_paired(obj, clt1, clt2, paste0("6h_", as.character(clt1), "_",  as.character(clt2)))
```

## Pyscenic expression network inference

Network inference with machine learning

```python
#########################################################
# PART1
# arboreto and GRNBoost2
#########################################################
# https://academic.oup.com/bioinformatics/article/35/12/2159/5184284
import pandas as pd
import numpy as np

from arboreto.algo import grnboost2
from distributed import LocalCluster, Client

import sys

if __name__=='__main__':

    local_cluster = LocalCluster(n_workers=10,
                                threads_per_worker=1,
                                memory_limit=8e9)
    custom_client = Client(local_cluster)
    
    ex_matrix = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=-1)
    ex_matrix.shape
    tf_df = pd.read_csv(sys.argv[2], sep="\t", header=0, index_col=None)
    tf_df['tf'] = tf_df.apply(lambda x: x['gene'].replace("A157_", ""), axis=1)
    tf_list = list(tf_df['tf'])
    
    adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_list, verbose=True, seed=12345)
    adjancencies.to_csv(sys.argv[3], index=False, sep='\t')
```

TF binding motif identification

```python
#########################################################
# PART2
# cisTarget
#########################################################

#########################################################
# STEP1 motif meme -> cluster burst
# python
#########################################################

from Bio import motifs

motifs_list = []
motifs_f = open(r"Stu_TF_binding_motifs.meme", "r")
for m in motifs.parse(motifs_f, "minimal"):
    motifs_list.append(m)
print(motifs_list[0].counts)
print(dir(motifs_list[0]))
motifs_f.close()

gene_dict = {}
with open(r"A157mapDM.txt", "r") as map_f:
    for line in map_f:
        tmp = line.strip("\n").split("\t")
        potato_g = tmp[0].split('.')[0]
        DM_g = tmp[1]
        if DM_g in gene_dict.keys():
            if potato_g not in gene_dict[DM_g]:
                gene_dict[DM_g].append(potato_g)
            else:
                pass
        else:
            gene_dict[DM_g] = [potato_g]
            
for m in motifs_list:
    tmp_m = m
    if m.name not in gene_dict.keys():
        continue
    for g in gene_dict[m.name]:
        with open(g.replace("A157_", "") + ".cb", "w") as outf:
            tmp_m.name = g
            motif_cb = motifs.write([tmp_m], "clusterbuster")
            outf.write(motif_cb)
```

```shell
#########################################################
# STEP2 prepare cisTarget db
# shell
#########################################################
ls Stu_TF_binding_motifs_cb | sed 's/\.cb//g' > motifs.lst

python /path/to/create_cisTarget_databases/create_cisTarget_databases-master/create_cistarget_motif_databases.py \
-f A157_franking_10kb.rename.fa \
-M Stu_TF_binding_motifs_cb \
-m motifs.lst \
-o A157 \
-t 20 \
-c cbust \
-g #[0-9]+$ # 必须 通过motifs id 获取gene id 
```

```python
#########################################################
# STEP3 motif ann
# python
#########################################################

import os
import pandas as pd

flist = [f.replace(".cb", "") for f in os.listdir("Stu_TF_binding_motifs_cb")]
df = pd.DataFrame(
    {"#motif_id": flist}
)
df["gene_name"] = df["#motif_id"]
df["motif_similarity_qvalue"] = 0
df["orthologous_identity"] = 1
df["description"] = df["#motif_id"]

df.to_csv("motif.tbl", header=True, index=False, sep='\t')
```

Pyscenic network inference

```python
#########################################################
# PART3
# run pyscenic
#########################################################

#########################################################
# STEP1 prepare_ctx_input.py
# python
#########################################################
import pandas as pd

# exp mat
EXP_MAT_FNAME = "exp_mat.tab"

# RIBOSOMAL genes
ribo_df = pd.read_csv("Ribo.txt", header=0, index_col=None, sep="\t")
RIBOSOMAL_GENE = list(ribo_df['gene'])

# adjancencies matrix from GRNBoost2
ADJACENCIES_FNAME = "adjancencies.tsv"
ADJACENCIES_FILTED_FNAME = "adjancencies_filted.tsv"

# motif from cisTarget
DATABASES_GLOB = "A157.genes_vs_motifs.rankings.feather"
MOTIF_ANNOTATIONS_FNAME="motif.tbl"

def edges_filter(adj_df):
    # Derive potential regulomes from these co-expression modules 
    adj_sorted_df = adj_df.sort_values(['TF', 'importance'], ascending=[1,0])

    # remove ribosomal
    adj_sorted_filted_df = adj_sorted_df[(~adj_sorted_df['target'].isin(RIBOSOMAL_GENE))]

    return adj_sorted_filted_df
    
# RIBOSOMAL
RIBOSOMAL_GENE = []

if __name__=='__main__':

    ex_matrix = pd.read_csv(EXP_MAT_FNAME, sep='\t', header=0, index_col=-1)
    ex_matrix.to_csv("exp_mat.tsv", sep="\t", header=True, index=True) # cell x genes 3000 x36000

    adjacencies = pd.read_csv(ADJACENCIES_FNAME, header=0, index_col=None, sep='\t')
    assert adjacencies.shape[1] == 3, "check adjacenies whether it is right"
    adj_sorted_filted_df = edges_filter(adjacencies)
    adj_sorted_filted_df.columns = ["TF", "target", "importance"]
    adj_sorted_filted_df.to_csv(ADJACENCIES_FILTED_FNAME, header=True, index=False, sep="\t")
```

```shell
#########################################################
# STEP2 run_ctx
# shell
#########################################################

# 注意输入
# 邻接矩阵，表达矩阵必须以csv或tsv结尾
/path/to/pyscenic ctx adjancencies.tsv A157.genes_vs_motifs.rankings.feather \
--annotations_fname motif.tbl \
--expression_mtx_fname ./exp_mat.tsv \
--top_n_targets 50 \ # 筛选邻接边的方法
--num_workers 20 \
--nes_threshold 2 \ # module聚类阈值
--mode "dask_multiprocessing" \
--output reg.csv \
--no_pruning \ # 不进行过滤
--min_genes 3

```



