
scRNA <- scRNA_celltype
table(scRNA$majorClass)
cell_selection <- subset(scRNA, subset = majorClass=="Astrocytes") 
table(cell_selection$group)

options(stringsAsFactors = F)
options(future.globals.maxSize = 20000 * 1024^2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(data.table)
library(harmony)
library(clustree)

clean_cell <- cell_selection@meta.data
clean_count <- table(clean_cell$orig.ident)
#C2 C3-1 C3-2   L1   L2   L3 
#1388  650 1144 1352 1811 1329 
clean_count
sum(clean_count)
#[1] 2439
scRNA_harmony <- NormalizeData(cell_selection) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose=FALSE)
plot1 <- VariableFeaturePlot(scRNA_harmony)
plot2 <- LabelPoints(plot=plot1, points=head(VariableFeatures(scRNA_harmony),10), repel=TRUE, size=2.5)
plot2
ggsave(plot2, file="VariableFeatures.pdf", width=6, height=6)

DimPlot(scRNA_harmony, reduction="pca", group.by="orig.ident")
ggsave( file="dimplot.pdf", width=6, height=6)
elbowplot <- ElbowPlot(scRNA_harmony, ndims=50)  #30
elbowplot
ggsave(elbowplot, file="elbowplot.pdf", width=6, height=6)

pc.nums=13
scRNA_harmony <- JackStraw(scRNA_harmony, num.replicate=100, dims=pc.nums)
scRNA_harmony <- ScoreJackStraw(scRNA_harmony, dims=1:pc.nums)
jackstrawplot <- JackStrawPlot(scRNA_harmony, dims=1:pc.nums)
jackstrawplot
ggsave(jackstrawplot, file="jackstrawplot.pdf", width=6, height=6)

system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars="orig.ident")})

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction="harmony", dims=1:pc.nums) 
scRNA_harmony <- FindClusters(
  object=scRNA_harmony,
  resolution=c(seq(.1,1,.1))
)  
clustreeplot <- clustree(scRNA_harmony@meta.data, prefix="RNA_snn_res.")  #0.4
clustreeplot
ggsave(clustreeplot, file="clustreeplot.pdf", width=12, height=10)


scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction="harmony", dims=1:pc.nums)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution=0.4) 
table(scRNA_harmony@active.ident)  
#0    1    2    3    4    5 
#3427 2682 1257  498  181  111 
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction="harmony", dims=1:pc.nums)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction="harmony", dims=1:pc.nums)

plot1 <- DimPlot(scRNA_harmony, reduction="umap", group.by="orig.ident")
#plot3 <- DimPlot(scRNA_harmony, reduction="tsne", label=T)
plot2 <- DimPlot(scRNA_harmony, reduction="umap", label=T)
plot1
plot2
ggsave(plot1, file="umapplot1.pdf", width=8, height=6)
ggsave(plot2, file="umapplot2.pdf", width=8, height=6)
saveRDS(scRNA_harmony, "scRNA_AS.rds")

scRNA_harmony <- macro
scRNA_harmony <- readRDS("scRNA_AS.rds")
DotPlot(scRNA_harmony, features=c("AQP4", "SLC1A2")) + RotatedAxis()  #Astrocytes,1，5
table(scRNA_harmony$seurat_clusters)
scRNA <-subset(scRNA_harmony,seurat_clusters!=8)
table(scRNA$seurat_clusters)
scRNA_harmony <- scRNA

table(scRNA$seurat_clusters)
DotPlot(scRNA_harmony, features=c("CDKN2A","CDKN1A","CHI3L1","THBS2","HGF","S100B","MALT1","ENO1")) ##mono  2
DotPlot(scRNA_harmony, features=c("Adgre1", "Cd68", "Mrc1")) ##Macrophages 0,1,3,4,5,7
DotPlot(scRNA_harmony, features=c("Itgax", "H2-Ab1")) ##Dendritic Cells 6
DotPlot(scRNA_harmony, features=c("Ly6g", "Cd47", "S100a8", "S100a9")) ##Neutrophils

scRNA_harmony@meta.data$minor <- NA
scRNA_harmony@meta.data$minor <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(0),"Astrocytes1",scRNA_harmony@meta.data$minor)
scRNA_harmony@meta.data$minor <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(1),"Astrocytes2",scRNA_harmony@meta.data$minor)
scRNA_harmony@meta.data$minor <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(2),"Astrocytes3",scRNA_harmony@meta.data$minor)
scRNA_harmony@meta.data$minor <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(3),"Astrocytes4",scRNA_harmony@meta.data$minor)
scRNA_harmony@meta.data$minor <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(4),"Astrocytes5",scRNA_harmony@meta.data$minor)
scRNA_harmony@meta.data$minor <- ifelse(scRNA_harmony@meta.data$seurat_clusters%in%c(5),"Astrocytes6",scRNA_harmony@meta.data$minor)
scRNA_celltype <- scRNA_harmony
table(scRNA_harmony$minor)
Idents(scRNA_harmony) <- "minor"
DimPlot(scRNA_harmony,group.by = "group")
DimPlot(scRNA_harmony,group.by = "minor")


Idents(scRNA_celltype) <- "minor"  

degs2 <- FindAllMarkers(scRNA_celltype, logfc.threshold=0.5,
                        test.use="wilcox", 
                        return.thresh=0.25, 
                        min.pct=0.3, only.pos=F)
write.csv(degs2, file="marker_wilcox.csv")


DimPlot(scRNA_harmony,group.by = "minor",label = T)

ggsave( file="umapplot3.pdf", width=8, height=6)

cells <- c("Monocytes","Macrophages","Dendritic Cells")
genes_to_check <- c("Ccr2","Cd14",
                    "Adgre1", "Cd68", "Mrc1",
                    "Itgax", "H2-Ab1"
)
scRNA_celltype <- scRNA_harmony
scRNA_celltype$minor <- factor(scRNA_celltype$minor,levels=cells)
DotPlot(scRNA_celltype,features=genes_to_check,
        assay='RNA',group.by = "minor")+coord_flip()+ 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL)+ 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
ggsave("bubble.pdf",width=6, height=4)


cells<- c("Astrocytes1","Astrocytes2","Astrocytes3","Astrocytes4","Astrocytes5","Astrocytes6")
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
          legend.position=c(0.5,0.9),
          topGeneN = 5)
ggsave("Volcan.pdf", width=10, height=6)

library(Seurat)
library(ClusterGVis)
library(org.Hs.eg.db)
library(COSG)
library(dior)
library(dplyr)
degs2$cluster <- factor(degs2$cluster, levels=cells)
sce.markers <- degs2 %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n=30, wt=avg_log2FC)


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


scRNA_celltype$minor <- factor(scRNA_celltype$minor,levels=cells)
#scRNA_celltype$majorClass <- factor(scRNA_celltype$majorClass,levels=cells)
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
           cluster.order = c(1:6) )
dev.off()


qual_col_pals <- brewer.pal.info[brewer.pal.info$category=="qual",]  
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pB_df <- table(scRNA_celltype@meta.data$minor, scRNA_celltype@meta.data$orig.ident) %>% reshape2::melt()
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

pB_df <- table(scRNA_celltype@meta.data$minor, scRNA_celltype@meta.data$group) %>% reshape2::melt()
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

pB1 <- pB1 + 
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
  theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.title.x = element_blank()) # 

pB2 <- pB2 + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
  theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.title.x = element_blank()) # 

# 调整图像的组合
Feature_ber <- ggarrange(pB1, pB2, ncol = 2, nrow = 1, heights = c(1, 1), 
                         common.legend = TRUE, legend = "right")
Feature_ber

ggsave("vis2.pdf", Feature_ber, width=12, height=6)







