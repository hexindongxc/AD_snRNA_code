
options(stringsAsFactors = F)
options(future.globals.maxSize = 20000 * 1024^2)

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(scales)

p <- DimPlot(scRNA_harmony, reduction="umap", group.by="seurat_clusters", label=T) 
p
DotPlot(scRNA_harmony, features="PTPRC") + RotatedAxis()  #non-immune,0,1,5,6,7,8,11
DotPlot(scRNA_harmony, features=c("EPCAM","KRT18","KRT19","CD24"))  #epithelial,2,3,4,12
DotPlot(scRNA_harmony, features=c("MME","COL1A1","DCN","FGF7"))  #fibo,9,13
DotPlot(scRNA_harmony, features=c("CCL21","PDPN"))  #fibroblastic reticular cell
DotPlot(scRNA_harmony, features=c("CLDN5","CDH5","PECAM1","VWF")) + RotatedAxis()  #endo,10
DotPlot(scRNA_harmony, features=c("RGS5","PDGFRB","ACTA2")) + RotatedAxis()  #Mural cells,
DotPlot(scRNA_harmony, features=c("RGS5","MCAM","ACTA2")) + RotatedAxis()  #pericyte
DotPlot(scRNA_harmony, features=c("CD44","PROM1","LGR5")) + RotatedAxis()  #CSCs
DotPlot(scRNA_harmony, features=c("COL1A1","COL3A1")) + RotatedAxis()  #mesenchymal stromal cells(MSCs)，
DotPlot(scRNA_harmony, features=c("DES","SMOC2","KCNMA1")) + RotatedAxis()    #myoepithelial cells
DotPlot(scRNA_harmony, features=c("NRXN1","CDH19","PCSK2")) + RotatedAxis()    #Neuron like cell,PMID: 33397933
DotPlot(scRNA_harmony, features=c("PLN","HIGD1B","ITGA7")) + RotatedAxis()    #smooth muscle cell，

DotPlot(scRNA_harmony, features=c("HLA-DRA", "CX3CR1", "C1QB", "CSF1R","CD74")) + RotatedAxis()  #Microglia 4
DotPlot(scRNA_harmony, features=c("AQP4", "SLC1A2")) + RotatedAxis()  #Astrocytes,1，5
DotPlot(scRNA_harmony, features=c("SYT1", "GRIK2", "GRIA1", "GRIN2B", "RBFOX1")) + RotatedAxis()  #Neurons 3,7,9 
DotPlot(scRNA_harmony, features=c("PCDH15", "MEGF11")) + RotatedAxis()  #OPCs,2
DotPlot(scRNA_harmony, features=c("FLT1", "CLDN5")) + RotatedAxis()  #Endothelial Cells，8
DotPlot(scRNA_harmony, features=c("MOBP", "MBP", "PLP1")) + RotatedAxis()  #Oligodendrocytes,0,6
DotPlot(scRNA_harmony, features=c("CDH18","CADPS2")) + RotatedAxis()  #
table(scRNA_harmony$seurat_clusters)
Idents(scRNA_harmony) <- "seurat_clusters"  
DimPlot(scRNA_harmony)
degs2 <- FindAllMarkers(scRNA_harmony, logfc.threshold=0.5,
                        test.use="wilcox", 
                        return.thresh=0.25, 
                        min.pct=0.3, only.pos=F)
write.csv(degs2, file="marker_wilcox1.csv")

scRNA_harmony@meta.data$majorClass <- NA
scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(4),"Microglia",scRNA_harmony@meta.data$majorClass)
scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(1,5),"Astrocytes",scRNA_harmony@meta.data$majorClass)
scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(3,7,9),"Neurons",scRNA_harmony@meta.data$majorClass)
#scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters==/,"Mural cells",scRNA_harmony@meta.data$majorClass)
#scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters==/,"Neuron like cells",scRNA_harmony@meta.data$majorClass)

scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(2),"OPCs",scRNA_harmony@meta.data$majorClass)
#scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(6),"NK cells",scRNA_harmony@meta.data$majorClass)
#scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(0),"CD4+T cells",scRNA_harmony@meta.data$majorClass)
scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(8),"Endothelial Cells",scRNA_harmony@meta.data$majorClass)
scRNA_harmony@meta.data$majorClass <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(0,6),"Oligodendrocytes",scRNA_harmony@meta.data$majorClass)

p2 <- DimPlot(scRNA_harmony, group.by="majorClass", label=T, label.size=5, reduction="umap", repel=TRUE)
p2
#ggsave("celltypeumap_tsne.pdf", p1, width=6, height=6)
ggsave("celltypeumap_umap.pdf", p2, width=8, height=6)
saveRDS(scRNA_harmony, "scRNA_celltype.rds")

#######################################################################################

scRNA_celltype<-scRNA_harmony
scRNA_celltype$celltype <- scRNA_celltype$majorClass
cells <- c("Oligodendrocytes","Astrocytes","OPCs","Microglia","Neurons","Endothelial Cells")
genes_to_check <- c("PLP1","MBP","MOBP",  
                    "AQP4", "SLC1A2", 
                    "PCDH15", "MEGF11",
                    "CD74","HLA-DRA", "CX3CR1", 
                    "SYT1", "GRIK2", "RBFOX1",
                    "DCN","FLT1", "CLDN5"
)
scRNA_celltype$majorClass <- factor(scRNA_celltype$majorClass,levels=cells)

DotPlot(scRNA_celltype,features=genes_to_check,
        assay='RNA',group.by = "majorClass")+coord_flip()+ 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL)+ 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 

ggsave("bubble.pdf",width=8, height=6)
library(scCustomize)
Idents(scRNA_celltype) <- "majorClass"
Clustered_DotPlot(seurat_object = scRNA_celltype, features = genes_to_check,cluster_ident=F,cluster_feature=F, exp_color_min = -1, exp_color_max=2,x_lab_rotate=90)


Idents(scRNA_celltype) <- "minor" 
scRNA_celltype$celltype <- factor(scRNA_celltype$celltype,levels=cells)
DimPlot(scRNA_celltype)

degs2 <- FindAllMarkers(scRNA_celltype, logfc.threshold=0.5,
                       test.use="wilcox", 
                       return.thresh=0.25, 
                       min.pct=0.3, only.pos=F)
write.csv(degs2, file="marker_minor_wilcox.csv")

Idents(scRNA_celltype) <- "group"  
DimPlot(scRNA_celltype)
degs2 <- FindAllMarkers(scRNA_celltype, logfc.threshold=0.5,
                        test.use="wilcox", 
                        return.thresh=0.25, 
                        min.pct=0.3, only.pos=F)
write.csv(degs2, file="marker_group_wilcox.csv")
library(scRNAtoolVis)
degs2$cluster <- factor(degs2$cluster, levels=cells)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
jjVolcano(diffData = degs2,
          log2FC.cutoff = 0.5,
          # col.type = "adjustP",
          tile.col = colour,
          size  = 3.5,
          celltypeSize=5,
          legend.position=c(0.75,0.9),
          topGeneN = 5)
ggsave("Volcan.pdf", width=12, height=8)

library(Seurat)
library(ClusterGVis)
library(org.Hs.eg.db)
library(COSG)
library(dior)
library(dplyr)
degs2$cluster <- factor(degs2$cluster, levels=cells)
sce.markers <- degs2 %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n=50, wt=avg_log2FC)


st.data <- prepareDataFromscRNA(object=scRNA_celltype,
                                diffData=sce.markers,
                                showAverage=TRUE)  
library(org.Hs.eg.db)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                      type = "BP",
                       organism = "hsa",
                      pvalueCutoff = 0.5,
                      topn = 6,
                     seed = 5201314)

sce.markers <- degs2 %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n=10, wt=avg_log2FC)

scRNA_celltype$celltype <- factor(scRNA_celltype$celltype,levels=cells)
scRNA_celltype$majorClass <- factor(scRNA_celltype$majorClass,levels=cells)
st.data <- prepareDataFromscRNA(object=scRNA_celltype,
                                diffData=sce.markers,
                                showAverage=TRUE)  

markGenes = sce.markers$gene
pdf('10DEGGpheatmap.pdf',height = 12,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           markGenes.side = "left",
           annoTerm.data = enrich,
           show_row_dend = F,
           markGenes = markGenes,
           cluster.order = c(1:6),
           boxcol =colour[0:100] )
dev.off()
qual_col_pals <- brewer.pal.info[brewer.pal.info$category=="qual",]  
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pB_df <- table(scRNA_celltype@meta.data$majorClass, scRNA_celltype@meta.data$orig.ident) %>% reshape2::melt()
colnames(pB_df) <- c("Cluster","Group","Number")
pB_df$Cluster <- factor(pB_df$Cluster,levels=cells)
sample_color <- col_vector[1:12]

pB1 <- ggplot(data=pB_df, aes(x=Group, y=Number, fill=Cluster)) +
  geom_bar(stat="identity", width=0.8, position="fill") +
  scale_fill_manual(values=col_vector[1:13]) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="", y="Ratio") +
  theme(axis.text.y=element_text(size=12, colour="black")) +
  theme(axis.text.x=element_text(size=12, colour="black")) +
  theme(axis.text.x.bottom=element_text(hjust=1, vjust=1, angle=45))
pB1

pB_df <- table(scRNA_celltype@meta.data$majorClass, scRNA_celltype@meta.data$group) %>% reshape2::melt()
colnames(pB_df) <- c("Cluster","Group","Number")
pB_df$Cluster <- factor(pB_df$Cluster,levels=cells)
pB_df$Group <- factor(pB_df$Group,levels=c("CT","AD"))
sample_color <- col_vector[1:12]

pB2 <- ggplot(data=pB_df, aes(x=Group, y=Number, fill=Cluster)) +
  geom_bar(stat="identity", width=0.8, position="fill") +
  scale_fill_manual(values=col_vector[1:13]) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="", y="Ratio") +
  theme(axis.text.y=element_text(size=12, colour="black")) +
  theme(axis.text.x=element_text(size=12, colour="black")) +
  theme(axis.text.x.bottom=element_text(hjust=1, vjust=1, angle=45))
pB2
pB1 <- pB1 + theme(legend.position = "none", plot.margin = unit(c(1, 1, 1, 1), "cm"))
pB2 <- pB2 + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
Feature_ber <- ggarrange(pB1, pB2, ncol = 2, nrow = 1, heights = c(1, 1))
print(Feature_ber)
ggsave("Vis1.pdf", Feature_ber, width=12, height=6)


library(sscVis)
library(data.table)
library(Startrac)
library(readr)
library(dplyr)
library(sscClust)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(patchwork)
library(stringr)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(scales)
library(ggalluvial)
library(ggpubr)
library(reshape2)
library(dior)

out.prefix = "./Roe"
if(T){
  text.size = 8
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   # axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(), 
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
                   # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
                   # strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
  )
}

scRNA_celltype@meta.data$NICB <- scRNA_celltype$group
cellInfo.tb <- scRNA_celltype@meta.data
cellInfo.tb$macro <- as.character(cellInfo.tb$macro)

startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient=F,
                                                 colname.cluster="majorClass",
                                                 colname.patient="orig.ident",
                                                 colname.tissue="NICB"))
Roe.result <- cbind(id=rownames(startrac.dist),startrac.dist)
cuts <- c(0.0000000, 0.0000001, 0.7,1 ,1.4,  2)
startrac.dist.bin.values <- factor(c("-", "+/-", "+", "++", "+++"),levels=c("-", "+/-", "+", "++", "+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                            ncol=ncol(startrac.dist))
colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist)
startrac.dist.bin

sscVis::plotMatrix.simple(startrac.dist,
                          out.prefix=sprintf("%s.startrac.dist",out.prefix),
                          show.number=F,
                          clust.row=T,
                          #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                          exp.name=expression(italic(R)[o/e]),
                          z.hi=2,
                          #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                          palatte=viridis::viridis(7),
                          pdf.width = 4.5, pdf.height = 4.5)

startrac.dist1 = ifelse(startrac.dist>2,2,startrac.dist)

plot_2 = pheatmap::pheatmap(startrac.dist1,
                            color = viridis::viridis(7),
                            cluster_rows = T,
                            cluster_cols = F,
                            fontsize = 10,
                            border_color = "white",
                            treeheight_col = 0,
                            treeheight_row = 0
)

plot_data = reshape2::melt(startrac.dist1)
head(plot_data)
colnames(plot_data) = c("meta.cluster", "loc", "Roe")
startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
# plot_data$bin = startrac.dist.bin2$value %>% factor(levels=c("-", "+/-", "+", "++", "+++"))
plot_data$bin = startrac.dist.bin2$value %>% factor(levels=c("-","+/-", "+", "++", "+++"))
plot_data$meta.cluster = factor(plot_data$meta.cluster,
                                levels = rev(plot_2$tree_row$labels[plot_2$tree_row$order]))
head(plot_data)
plot_data$meta.cluster <- factor(plot_data$meta.cluster,levels=cells)
plot_data$loc <- factor(plot_data$loc,levels=(c("CT","AD")))
p3 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = Roe)) +
  geom_tile() +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  scale_fill_gradientn(name = "Ro/e",
                       colours = viridis::viridis(7)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme
p3

p4 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = bin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Ro/e",
    values = c(
      "NA"="white",
      "-" = "white",
      "+/-" = "#fde6ce",
      "+" = "#fcc08b",
      "++" = "#f5904a",
      "+++" = "#e6540d"
    ),
    labels = c(0, 0.7,1 ,1.4, ">2")
  ) +
  geom_text(
    aes(label = bin),
    color = "black",
    size = 3,
    show.legend = T
  ) +
  guides(fill = guide_legend(title = "Ro/e",
                             override.aes = list(label = c(
                               
                               "+/-", "+","++","+++"
                             )))) + theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm")) + mytheme;p4

wrap_plots(p3,p4)
ggsave(filename = paste0(out.prefix,".P3_P4.ggplot.heatmap.pdf"), 
       height = 5, width = 7)

p4
ggsave(filename = "P4.ggplot.heatmap.pdf", 
       height = 5, width = 3.5)
