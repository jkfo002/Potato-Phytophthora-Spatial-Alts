# Giotto pipeline

```R
library(Giotto)
library(Seurat)
library(ggplot2)
set.seed(12345)
library(ComplexHeatmap)
#library(viridis)

args <- commandArgs(T)

dat_part=readRDS("final_used_data.merged.RDS")
dat_part=subset(dat_part, subset=hours==args[1])
#dat_part=subset(dat_part, subset=(hours=="24T"))
#dat_part <- FindVariableFeatures(dat_part, nfeautures=3000)

dir.create(args[1])

coor.dat <- dat_part@meta.data[,c("x", "y", "cellid")]
colnames(coor.dat) <- c("sdimx", "sdimy", "cell_ID")
meta.data <- dat_part@meta.data
meta.data$cell_ID <- meta.data$cellid

instrs = createGiottoInstructions(python_path = "/path/to/envs/giotto/bin/python3",
                                  save_dir = getwd(),
                                  plot_format = 'png',
                                  dpi = 200,
                                  height = 9,
                                  width = 9)

gobject = createGiottoObject(
    expression = list(
        raw=dat_part@assays$Spatial@counts, 
        normalized=dat_part@assays$Spatial@data,
        scaled=dat_part@assays$Spatial@scale.data
    ), 
    spatial_locs = coor.dat, # neceessary
    instructions = instrs
    )
gobject <- addCellMetadata(
    gobject, new_metadata=meta.data,
    by_column = T, column_cell_ID = "cellid"
)

# test plot
#p <- spatPlot(gobject,
#         cell_color = 'seurat_clusters',
#         return_plot = TRUE)
#ggsave(plot=p, filename="./plot/test.pdf", width=50, height=50, limitsize=F)  

# network
# calculate spatial variant genes same with SpatialIDE
## create spatial network
gobject <- createSpatialNetwork(gobject)
p <- spatPlot2D(gobject = gobject,
           show_network= T,
           point_size=1,
           network_color = 'blue',
           return_plot = TRUE) # ggplot2 obj
ggsave(plot=p, filename=paste0("cell_network.", args[1], ".pdf"), width=20, height=20, limitsize=F)
## test spatial variation network    
ranktest = binSpect(gobject, bin_method = 'rank',
                    calc_hub = T, hub_min_int = 5,
                    spatial_network_name = 'Delaunay_network')
coexpression.gene <- VariableFeatures(dat_part)
input_gene <- ranktest[1:2000,]$feats

## calculate pairwise distances between genes (but set network_smoothing=0 to use default clustering)
spat_cor_netw_DT = detectSpatialCorFeats(gobject,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network',
                                         subset_feats = input_gene)
# spatial related genes failed
# top10_genes = showSpatialCorFeats(spat_cor_netw_DT, feats = '11G014750', show_top_feats = 10)

## cluster
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, name = 'spat_netw_clus', k = 10)
p <- heatmSpatialCorFeats(gobject,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus',
                     return_plot = T,
                     col = colorRamp2(seq(0,0.3, length=9), rev(brewer.pal(9, "RdBu")))) # 可以加上ComplexHeatmap的参数
pdf(paste0("heatmap.", args[1], ".k10.pdf"), width = 20, height = 20)
p <- draw(p)
dev.off()

clt.tab <- as.data.frame(spat_cor_netw_DT$cor_clusters)
colnames(clt.tab) <- c("modules")
clt.tab$gene <- rownames(clt.tab)
write.table(clt.tab, file = paste0("modules.", args[1], ".tab"), col.names = T, row.names = F, sep = "\t", quote = F)

save.image(file=paste0("giotto.", args[1], ".RData"))

#######################################################
## iterated find spatial neighbor
#######################################################
library(Giotto)
library(Seurat)
library(ggplot2)
set.seed(12345)
library(ComplexHeatmap)
library(viridis)

my_color <- c("l1"="#a3cd9e", "l2"="#529471", "l3"="#35635b", "pi"="#91684a", "others"="grey")

load("inf.giotto.20231020.RData")
pi_cell <- read.table("inf_pi.tab", header=T, sep="\t")
pi_cell <- pi_cell[pi_cell$nfeature_10.12=="Y",][["cellid"]]

pi.neigbor <- findNetworkNeighbors(gobject, spatial_network_name="Delaunay_network", source_cell_ids = pi_cell)
pi_cells_l1 <- pi.neigbor[pi.neigbor$nb_cells=="neighbor", ][["cell_ID"]]
pi.neigbor_s2 <- findNetworkNeighbors(gobject, spatial_network_name="Delaunay_network", source_cell_ids = c(pi_cell, pi_cells_l1))
pi_cells_l2 <- pi.neigbor_s2[pi.neigbor_s2$nb_cells=="neighbor", ][["cell_ID"]]
pi.neigbor_s3 <- findNetworkNeighbors(gobject, spatial_network_name="Delaunay_network", source_cell_ids = c(pi_cell, pi_cells_l1, pi_cells_l2))
pi_cells_l3 <- pi.neigbor_s3[pi.neigbor_s3$nb_cells=="neighbor", ][["cell_ID"]]

out <- as.data.frame(pi.neigbor)
out$pi_layers <- "others"
out$cellid <- out$cell_ID
out$pi_layers <- ifelse(out$cellid %in% pi_cell, "pi", out$pi_layers)
out$pi_layers <- ifelse(out$cellid %in% pi_cells_l1, "l1", out$pi_layers)
out$pi_layers <- ifelse(out$cellid %in% pi_cells_l2, "l2", out$pi_layers)
out$pi_layers <- ifelse(out$cellid %in% pi_cells_l3, "l3", out$pi_layers)

potato.merged <- readRDS("final_used_data.merged.RDS")
potato.merged <- subset(potato.merged, subset=inf_stat=="Inf")
potato.merged <- subset(potato.merged, subset=cellid %in% out$cellid)

dat.plot <- as.data.frame(potato.merged@meta.data)
dat.plot <- left_join(dat.plot, out, by = "cellid")
draw_sample_combine=function(dat.plot, fpath) {
  #dat.plot <- tissue_dat[tissue_dat$gene_exp>0,]
  #dat.plot.0 <-  tissue_dat[tissue_dat$gene_exp==0,]
  p <- ggplot() +
    #  16是实心圆，21是空心圆, 21才可以画出边界
    geom_point(data=dat.plot,mapping=aes(x=y,y=x, color=pi_layers), size=0.5, shape=16) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_rect(fill="black"),plot.background=element_rect(fill="black")) +
    scale_color_manual(values = my_color) +
    #scale_color_gradientn(colours =c("#555B73","#00F2FF","#1CFF00","#FFF500","#FF4700")) +
    #scale_color_gradientn(fill = Greens) +
    #scale_color_gradientn(colours =c("#fafcc8", "#ff0000")) +
    #scale_fill_distiller(palette = "Greens") +
    #scale_color_gradientn(colours =c("#241333","#312547","#8B1C1D","#FAEA00")) +
    coord_fixed(ratio=1)

  ggsave(fpath, p, width=50, heigh=50, limitsize = FALSE)
  
}
draw_sample_combine(dat.plot, "./layers.pdf")
```

