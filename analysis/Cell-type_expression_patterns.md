# Cell-type expression patterns

## DEG

DEG identification

```R
library(tidyverse)
library(ggplot2)
library("dplyr")
library("Seurat")

# /public/home/daijichen/01_workdir/02_liyuying_2/00_fig4/fig_leaf_tissue/fig2/fig2/03_cell_type_heatmap/05_dotplot
#potato_merged=readRDS('/public/home/daijichen/01_workdir/02_liyuying_2/final_used_data_0616.merged.RDS')
#/home/dongzhaonian/anaconda3/envs/seurat_405/bin/R

my_color = c(
    "0"="#93d0fc", "1"="#2aaae1", "2"="#3f5246", "3"="#ad3074", "4"="#E7B10A", 
    "5"="#91a3b0", "6"="#d58282", "7"="#418557", "8"="#6dc831", "9"="#eacd01",
    "10"="#ff9f22", "11"="#af5b2b", "12"="#2e5073", "13"="#E6D2AA", "14"="#42032C",
    "15"="#752f9a", "16"="#721c28", "17"="#2c2cff", "18"="#d1bac8", "19"="#009580",
    "20"="#94ac78", "21"="#9f59f8", "22"="#ffff00", "23"="#fc0cd4", "24"="#ffccff"
)

# celltype changed
potato_merged=readRDS('/home/dongzhaonian/job/00_ST/06_potato_4/00_fig/marker_gene_dot/final_used_data_0807.merged.RDS')

# rds change
ct_rds_list <- list()

#haimian
dat1=subset(potato_merged,subset=leaf_tis=='haimian')
dat2=subset(dat1,idents=c(3,4,6,10,16,0,1,2,12))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["haimian"]] = dat3

#zhalan
dat1=subset(potato_merged,subset=leaf_tis=='zhalan')
dat2=subset(dat1,subset=(seurat_clusters==11 | (hours=='3T' & seurat_clusters=='7')|
                        (hours=='6h' & seurat_clusters=='7')|
                        (hours=='24T' & seurat_clusters=='7')|
                        (hours=='12T' & seurat_clusters=='17')))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["zhalan"]] = dat3

#shangbiaopi
dat1=subset(potato_merged,subset=leaf_tis=='shangbiaopi')
dat2=subset(dat1,idents=c(3,4,14,0,1,5))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["shangbiaopi"]] = dat3

#xiabiaopi
dat1=subset(potato_merged,subset=leaf_tis=='xiabiaopi')
dat2=subset(dat1,idents=c(6,14,16,0,5))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["xiabiaopi"]] = dat3

# weiguan
dat1=subset(potato_merged,subset=leaf_tis=='weiguan')
dat2=subset(dat1,idents=c(13,9,15))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["weiguan"]] = dat3

# baowei
dat1=subset(potato_merged,subset=leaf_tis=='baowei')
dat2=subset(dat1,idents=c(24,22))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["baowei"]] = dat3

diff_exp_paired <- function(obj, id1, id2, outprefix){
    diff.wilcox = FindMarkers(obj, only.pos = FALSE,
                                min.pct = 0.25, logfc.threshold = 0.1,
                                ident.1 = id1, ident.2 = id2
                                )
    diff.wilcox$gene = rownames(diff.wilcox)
    all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
    all.markers$gene=paste0('A157_',all.markers$gene)
    write.table(all.markers,
    file=paste0(outprefix, "_DEG_gene.txt"),
    sep="\t",col.names=T,quote=F,row.names=F)
}


for (n in names(obj_list)){
    dat_tmp <- obj_list[[n]]
    diff_exp_paired(dat_tmp, paste0("Inf",'_', n), paste0("CK",'_', n), paste0("Inf_vs_CK.", t, ".", n))
}
```

statistic

```R
library(reshape2)
library(ggVennDiagram)
library(ggplot2)
library(cowplot)
library(tidyverse)

d.list <- map(
  c("3T", "6h", "12T", "24T"), 
  function(x){
    dat <- read.table(paste0("celltype.DEG.", x, ".up_down_0.58.tab"), header = T, sep = "\t")
    colnames(dat) <- c("gene", "tissue", "time", "p")
    up_dat <- dat[dat$p=="up",]
    tmp <- as.data.frame(table(up_dat[["tissue"]]))
    colnames(tmp) <- c("cu", "count")
    tmp[["time"]] <- x
    
    return(tmp)
  }
)
up_plot.dat <- do.call(rbind, d.list)
up_plot.dat[["cu"]] <- factor(up_plot.dat[["cu"]], levels = c("shangbiaopi", "zhalan", "haimian", "xiabiaopi", "baowei", "weiguan"))
up_plot.dat[["time"]] <- factor(up_plot.dat[["time"]], levels = c("3T", "6h", "12T", "24T"))
up_plot.dat[["pattern"]] <- "up"

d.list <- map(
  c("3T", "6h", "12T", "24T"), 
  function(x){
    dat <- read.table(paste0("celltype.DEG.", x, ".up_down_0.58.tab"), header = T, sep = "\t")
    colnames(dat) <- c("gene", "tissue", "time", "p")
    up_dat <- dat[dat$p=="down",]
    tmp <- as.data.frame(table(up_dat[["tissue"]]))
    colnames(tmp) <- c("cu", "count")
    tmp[["time"]] <- x
    
    return(tmp)
  }
)
down.plot.dat <- do.call(rbind, d.list)
down.plot.dat[["cu"]] <- factor(down.plot.dat[["cu"]], levels = c("shangbiaopi", "zhalan", "haimian", "xiabiaopi", "baowei", "weiguan"))
down.plot.dat[["time"]] <- factor(down.plot.dat[["time"]], levels = c("3T", "6h", "12T", "24T"))
down.plot.dat[["pattern"]] <- "down"

plot.dat <- rbind(up_plot.dat, down.plot.dat)
color_manual <- c(
  "shangbiaopi"="#ffd166", "zhalan"="#06d6a0", "haimian"="#118ab2",
  "xiabiaopi"="#ef476f", "baowei" = "#caadff", "weiguan"="#f78c6b", "comm"="black"
)
top.mar=0.2
right.mar=0.1
bottom.mar=0.2
left.mar=0.1
mytheme<-theme_classic()+
theme(text=element_text(family = "sans",colour ="gray30",size = 12),
      axis.text.y = element_text(size = 12),
      axis.line.y=element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_line(size = 0.6,colour = "gray30"),
      axis.ticks = element_line(size = 0.6,colour = "gray30"),
      axis.ticks.length = unit(1.5,units = "mm"),
      plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                       units="inches"))

up_theme <- mytheme +  theme(axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())

p <- ggplot(up_plot.dat, 
            aes(x = time, 
                y = ifelse(pattern=="up",count,-count), 
                fill=cu, order=cu)) + 
  geom_bar(stat = "identity", position = 'dodge') +
  scale_fill_manual(values = color_manual) + up_theme  + scale_y_continuous(limits = c(0,510))
p1 <- ggplot(down.plot.dat, 
             aes(x = time, 
                 y = ifelse(pattern=="up",count,-count), 
                 fill=cu, order=cu)) + 
  geom_bar(stat = "identity", position = 'dodge') +
  scale_fill_manual(values = color_manual) + mytheme  + scale_y_continuous(limits = c(-510,0))
```

## Immune score

```R
# add mock
library(Seurat)
library(ggplot2)
library(tidyverse)
library(viridis)

# Specific gene list
gene_list <- read.table("gene.list", header = T, sep = "\t")[['gene']] %>% gsub("A157_", "", .)
# fill
my_color2 = c(
    "shangbiaopi"="#ffd166", "zhalan"="#06d6a0", "haimian"="#118ab2",
    "xiabiaopi"="#ef476f", "baowei" = "#caadff", "weiguan"="#f78c6b",
    "CK_shangbiaopi"="white", "CK_zhalan"="white", "CK_haimian"="white",
    "CK_xiabiaopi"="white", "CK_baowei" = "white", "CK_weiguan"="white"
)
# color
my_color3 = c(
    "CK_shangbiaopi"="#ffd166", "CK_zhalan"="#06d6a0", "CK_haimian"="#118ab2",
    "CK_xiabiaopi"="#ef476f", "CK_baowei" = "#caadff", "CK_weiguan"="#f78c6b",
    "shangbiaopi"="black", "zhalan"="black", "haimian"="black",
    "xiabiaopi"="black", "baowei" = "black", "weiguan"="black"
)

potato_merged=readRDS('final_used_data.RDS')
dat_part <- subset(
    potato_merged,
    subset=((leaf_tis2=="haimian")&(seurat_clusters %in% c(3,4,6,10,16))) |
    ((leaf_tis2=='zhalan')&(seurat_clusters==11)) |
    ((leaf_tis2=="shangbiaopi")&(seurat_clusters %in% c(3,4,14))) |
    ((leaf_tis2=="xiabiaopi")&(seurat_clusters %in% c(6,14,16))) |
    ((leaf_tis2=="weiguan")&(seurat_clusters==13)) |
    ((leaf_tis2=="baowei")&(inf_stat=="Inf")) | (inf_stat=="CK")
)
dat_part$leaf_tis3 <- ifelse(dat_part$inf_stat=="CK", paste0("CK_", dat_part$leaf_tis2), dat_part$leaf_tis2)
Idents(dat_part) <- dat_part$leaf_tis3
table(dat_part$leaf_tis3)

potato_merged <- AddModuleScore(potato_merged, features = list(gene_list), name="immune")
plot_df <- potato_merged@meta.data
plot_df$hours <- factor(plot_df$hours, levels = c("3T", "6h", "12T", "24T"))
plot_df$leaf_tis3 <- ifelse(plot_df$inf_stat=="CK", "CK", plot_df$leaf_tis2)
plot_df$leaf_tis3 <- factor(plot_df$leaf_tis3, levels = c("shangbiaopi", "CK_shangbiaopi", "zhalan", "CK_zhalan", "haimian", "CK_haimian", "xiabiaopi", "CK_xiabiaopi", "baowei", "CK_baowei", "weiguan", "CK_weiguan"))
plot_df <- plot_df[plot_df$cellid %in% dat_part$cellid,]

p <- ggplot(plot_df) +
geom_boxplot(aes(x=hours, y=immune1, fill=leaf_tis3, color=leaf_tis3)) +
scale_fill_manual(values = my_color2) +
scale_color_manual(values = my_color3) +
theme_bw() +
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave("celltype.immune_ms.pdf",p,width=12,height=5)

max_l <- c()
min_l <- c()
mean_l <- c()
group <- c()
time <- c()
for (i in c("shangbiaopi", "CK_shangbiaopi", "zhalan", "CK_zhalan", "haimian", "CK_haimian", "xiabiaopi", "CK_xiabiaopi", "baowei", "CK_baowei", "weiguan", "CK_weiguan")){
    for (j in c("3T", "6h", "12T", "24T")){
        tmp_df <- plot_df[plot_df[["leaf_tis3"]]==i & plot_df[["hours"]]==j,]
        max_l <- c(max_l, max(tmp_df[["immune1"]]))
        min_l <- c(min_l, min(tmp_df[["immune1"]]))
        mean_l <- c(mean_l, min(tmp_df[["immune1"]]))
        group <- c(group, i)
        time <- c(time, j)
    } 
}
plot_df2 <- data.frame(group=group, time=time, Max=max_l, Min=min_l, Mean=mean_l)
plot_df2[["group2"]] <- paste0(plot_df2[["group"]], "_", plot_df2[["time"]])
plot_df2[["group2"]] <- factor(
    plot_df2[["group2"]],
    levels = paste0(
        rep(c("shangbiaopi", "CK_shangbiaopi", "zhalan", "CK_zhalan", "haimian", "CK_haimian", "xiabiaopi", "CK_xiabiaopi", "baowei", "CK_baowei", "weiguan", "CK_weiguan"), times=4),
        "_",
        rep(c("3T", "6h", "12T", "24T"), each=12)
))

p2 <- ggplot(data=plot_df2,aes(x=group2,y=Mean))+
  geom_errorbar(aes(ymin=Min,ymax=Max),width=0, color=group)+
  geom_point(size=3, color=group)+
  scale_color_manual(values = my_color3)+
  theme_classic()+
  ylab("Module Score")
ggsave("celltype.immune_ms_line.20241006.pdf",p2,width=12,height=5)
```

## Immune pathway

```R
library(tidyverse)
library(ggplot2)
library("dplyr")
library("Seurat")

# /public/home/daijichen/01_workdir/02_liyuying_2/00_fig4/fig_leaf_tissue/fig2/fig2/03_cell_type_heatmap/05_dotplot
#potato_merged=readRDS('/public/home/daijichen/01_workdir/02_liyuying_2/final_used_data_0616.merged.RDS')
#/home/dongzhaonian/anaconda3/envs/seurat_405/bin/R

my_color = c(
    "0"="#93d0fc", "1"="#2aaae1", "2"="#3f5246", "3"="#ad3074", "4"="#E7B10A", 
    "5"="#91a3b0", "6"="#d58282", "7"="#418557", "8"="#6dc831", "9"="#eacd01",
    "10"="#ff9f22", "11"="#af5b2b", "12"="#2e5073", "13"="#E6D2AA", "14"="#42032C",
    "15"="#752f9a", "16"="#721c28", "17"="#2c2cff", "18"="#d1bac8", "19"="#009580",
    "20"="#94ac78", "21"="#9f59f8", "22"="#ffff00", "23"="#fc0cd4", "24"="#ffccff"
)

# celltype changed
potato_merged=readRDS('/home/dongzhaonian/job/00_ST/06_potato_4/00_fig/marker_gene_dot/final_used_data_0807.merged.RDS')

# rds change
ct_rds_list <- list()

#haimian
dat1=subset(potato_merged,subset=leaf_tis=='haimian')
dat2=subset(dat1,idents=c(3,4,6,10,16,0,1,2,12))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["haimian"]] = dat3

#zhalan
dat1=subset(potato_merged,subset=leaf_tis=='zhalan')
dat2=subset(dat1,subset=(seurat_clusters==11 | (hours=='3T' & seurat_clusters=='7')|
                        (hours=='6h' & seurat_clusters=='7')|
                        (hours=='24T' & seurat_clusters=='7')|
                        (hours=='12T' & seurat_clusters=='17')))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["zhalan"]] = dat3

#shangbiaopi
dat1=subset(potato_merged,subset=leaf_tis=='shangbiaopi')
dat2=subset(dat1,idents=c(3,4,14,0,1,5))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["shangbiaopi"]] = dat3

#xiabiaopi
dat1=subset(potato_merged,subset=leaf_tis=='xiabiaopi')
dat2=subset(dat1,idents=c(6,14,16,0,5))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["xiabiaopi"]] = dat3

# weiguan
dat1=subset(potato_merged,subset=leaf_tis=='weiguan')
dat2=subset(dat1,idents=c(13,9,15))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["weiguan"]] = dat3

# baowei
dat1=subset(potato_merged,subset=leaf_tis=='baowei')
dat2=subset(dat1,idents=c(24,22))
dat3=subset(dat2,subset=(inf_stat=='CK'|inf_stat=='Inf'))
dat3$new_group=paste0(dat3$inf_stat,'_',dat3$leaf_tis)
Idents(object = dat3)= dat3$new_group
ct_rds_list[["baowei"]] = dat3

# 对potato_merged通过各个ct的rds的cellid取子集

ct_cell_df <- data.frame(
    ct_cellid = c(
    rownames(ct_rds_list[["haimian"]]@meta.data), rownames(ct_rds_list[["zhalan"]]@meta.data), 
    rownames(ct_rds_list[["shangbiaopi"]]@meta.data), rownames(ct_rds_list[["xiabiaopi"]]@meta.data), 
    rownames(ct_rds_list[["weiguan"]]@meta.data)
    ),
    new_group = c(
        ct_rds_list[["haimian"]]@meta.data[['new_group']], ct_rds_list[["zhalan"]]@meta.data[['new_group']], 
        ct_rds_list[["shangbiaopi"]]@meta.data[['new_group']], ct_rds_list[["xiabiaopi"]]@meta.data[['new_group']], 
        ct_rds_list[["weiguan"]]@meta.data[['new_group']]
    )
)
dat.part.ct.violine <- subset(potato_merged, subset = cellid %in% ct_cell_df$ct_cellid)

# 对ct_cellid_set进行排序处理
ct_cell_df$ct_cellid <- factor(ct_cell_df$ct_cellid, levels = rownames(dat.part.ct.violine@meta.data))
ct_cell_df <- ct_cell_df[order(ct_cell_df$ct_cellid),]
dat.part.ct.violine$new_group <- ct_cell_df$new_group

# plot dot

Idents(dat.part.ct.violine) <- factor(dat.part.ct.violine$new_group, 
levels=c(paste0(rep(c("CK_", "Inf_"), each=5), rep(c("shangbiaopi", "zhalan", "haimian", "weiguan", "xiabiaopi"), time=2)))
)

plot_dotplot <- function(dat_part, gene.list, shortname, smax, ph, outputprefix){
    p1 <- DotPlot(object = dat_part, features = gene.list, scale.min = 0, scale.max=smax) + 
    theme_classic() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank()) +
    theme(axis.text.x = element_text(hjust = 1, vjust = .5,angle = 90,size=7),
    panel.background = element_rect(fill="white"),
    axis.title.x = element_text(colour = "white")) +
    coord_flip() + 
    scale_colour_gradient2(low = "lightgrey", high = "red") + 
    scale_x_discrete(labels = shortname) + scale_size(range=c(1,8))

    ggsave(p1, filename=outputprefix, width=6, height=ph)
}

gene.list <- read.table("./gene.list", header=T, sep="\t")

print(head(gene.list))

tmp.gene.list <- gene.list
feature_selected <- tmp.gene.list[['gene']] %>% unique(.) %>% gsub("A157_", "", .)
short_name <- tmp.gene.list[['rename']]

print(head(feature_selected))

plot_dotplot(
dat.part.ct.violine, 
feature_selected, short_name, 
50, 9,
paste0("./celltype_SAJA.20240129.pdf")
)
```

## Expression patttern

```R
set.seed(12345)
library("ggplot2")
library("dplyr")
library(tidyverse)
library(Matrix)
library(parallel)
library(reshape2)

my_color2 = c(
    "shangbiaopi"="#ffd166", "zhalan"="#06d6a0", "haimian"="#118ab2",
    "xiabiaopi"="#ef476f","weiguan"="#f78c6b"
)

# normalized expression tab
exp.dat <- read.table("exp.tab", header=T, sep="\t")
exp.melt.tab <- melt(exp.dat, id.vars = "gene")
colnames(exp.melt.tab) <- c("gene", "group", "value")
exp.melt.tab <- separate(exp.melt.tab, group, into = c("inf_stat", "celltype", "hours"), sep = "_")
exp.melt.tab$hours <- ifelse(exp.melt.tab$inf_stat=="CK", "0T", exp.melt.tab$hours)

reg <- "down"
genelist <- exp.melt.tab[['gene']] %>% unique(.)

for (n in names(my_color2)){
  if (!file.exists(paste0("./20230904_exp/", reg, "/", n))){
    dir.create(paste0("./20230904_exp/", reg, "/", n))
  }
}

gene_fc_lplot <- function(i){
  gene <- genelist[i]
  plot.tab <- exp.melt.tab[exp.melt.tab[['gene']]==gene,]
  plot.tab <- plot.tab[plot.tab[['celltype']]!='baowei',]
  if (reg == 'down'){
    plot.tab <- plot.tab[order(plot.tab[['value']], decreasing=F),]
  }else {
    plot.tab <- plot.tab[order(plot.tab[['value']], decreasing=T),]
  }
  
  ct_name <- plot.tab[['celltype']][1]

  plot.tab$hours <- factor(plot.tab$hours, levels=c("3T", "6h", "12T", "24T"))

  p <- ggplot(data = plot.tab, mapping = aes(x = hours, y = value, colour = celltype, shape = celltype, fill = celltype, group=celltype)) +
  geom_line(linetype="longdash", size=1) + geom_point(size=2.5) + theme_bw() + #绘制线图和点
  scale_color_manual(values = my_color2) +
  scale_shape_manual(values = c('zhalan'=7, 'weiguan'=9, 'haimian'=10, 'shangbiaopi'=12, 'xiabiaopi'=13))

  ggsave(plot = p, filename = paste0("./20230904_exp/", reg, "/", ct_name, "/", "A157_", gene, ".pdf"), 
  width = 7.5, height = 5)

}

mclapply(1:length(genelist), gene_fc_lplot, mc.cores=5)
```

