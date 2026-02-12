options(stringsAsFactors = F)
options(future.globals.maxSize = 20000 * 1024^2)

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
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(ggplot2)
library(ggpubr)

scRNA$majorClass <- factor(scRNA$majorClass,levels=cells)
macro <- scRNA_harmony
gene <- "SORBS1"

table(scRNA$majorClass)
table(macro$minor)
DimPlot(scRNA)

scRNA$majorClass <- factor(scRNA$majorClass,levels=cells)

macro <- scRNA_celltype

VlnPlot_scCustom(seurat_object = scRNA, features = gene,group.by = "majorClass",split.by = "group",colors_use = c("#1F77B4", "#FF7F0E"))+
  #stat_compare_means(method = "wilcox",label = "p.signif")+
 geom_boxplot(width=0.2,position=position_dodge(0.9), col="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("SORBS1_各类型细胞表达箱型图.pdf", width=10, height=5)

macro <-subset(macro,minor!="Astrocytes6")
macro$group <- factor(macro$group,levels=c("CT","AD"))
VlnPlot_scCustom(seurat_object = macro, features = gene,group.by = "minor",split.by = "group",colors_use = c("#1F77B4", "#FF7F0E"))+
  stat_compare_means(method = "wilcox",label = "p.signif")+
  geom_boxplot(width=0.2,position=position_dodge(0.9), col="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("SORBS1_各类型星状胶质细胞表达箱型图.pdf", width=10, height=5)


GeneGSEA <- function(scRNAdata, gene){
  GeneType <- as.data.frame(scRNAdata@assays$RNA$data[gene,])
  colnames(GeneType)[1] <- "genes"
  GeneType$group <- ifelse(GeneType$genes==0,"nonexpr","expr") 
  GeneType <- cbind(id=rownames(GeneType), GeneType)
  scRNAdata@meta.data$GeneGroup <- GeneType[match(rownames(scRNAdata@meta.data), GeneType$id), "group"]
  
  scRNAdata <- SetIdent(scRNAdata, value="GeneGroup")    
  DGE_cell_selection <- FindMarkers(scRNAdata, ident.2="nonexpr", ident.1="expr", min.pct=0.1, only.pos=FALSE)
  
  nrDEG <- DGE_cell_selection[,c("avg_log2FC","p_val")]
  colnames(nrDEG) <- c("log2FoldChange", "pvalue")
  
  gene <- bitr(rownames(nrDEG), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL, rownames(nrDEG))]
  
  geneList <- gene$logFC
  names(geneList) <- gene$ENTREZID
  geneList <- sort(geneList, decreasing=T)
  
  kk_gse <- gseKEGG(geneList=geneList, organism="hsa", minGSSize=1, pvalueCutoff=0.5, verbose=FALSE)
  kk_gse <- DOSE::setReadable(kk_gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
  sortkk <- kk_gse[order(kk_gse@result$enrichmentScore, decreasing=T),]
  
  results <- list(sortkk=sortkk, kk_gse=kk_gse)
  return(results)
}

GeneGSEA2 <- function(scRNAdata, gene){
  GeneType <- as.data.frame(scRNAdata@assays$RNA$data[gene,])
  colnames(GeneType)[1] <- "genes"
  GeneType$group <- ifelse(GeneType$genes> median(GeneType$genes),"expr","nonexpr")  
  GeneType <- cbind(id=rownames(GeneType), GeneType)
  scRNAdata@meta.data$GeneGroup <- GeneType[match(rownames(scRNAdata@meta.data), GeneType$id), "group"]
  
  scRNAdata <- SetIdent(scRNAdata, value="GeneGroup")    
  DGE_cell_selection <- FindMarkers(scRNAdata, ident.2="nonexpr", ident.1="expr", min.pct=0.1, only.pos=FALSE)
  
  nrDEG <- DGE_cell_selection[,c("avg_log2FC","p_val")]
  colnames(nrDEG) <- c("log2FoldChange", "pvalue")
  
  gene <- bitr(rownames(nrDEG), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL, rownames(nrDEG))]
  
  geneList <- gene$logFC
  names(geneList) <- gene$ENTREZID
  geneList <- sort(geneList, decreasing=T)
  
  kk_gse <- gseKEGG(geneList=geneList, organism="hsa", minGSSize=1, pvalueCutoff=0.5, verbose=FALSE)
  kk_gse <- DOSE::setReadable(kk_gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
  sortkk <- kk_gse[order(kk_gse@result$enrichmentScore, decreasing=T),]
  
  results <- list(sortkk=sortkk, kk_gse=kk_gse)
  return(results)
}
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(enrichR)
library(rafalib)
library(ggpubr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
TumorMacrophage <- macro

TumorMacrophage <- scRNA

GSEAResults <- GeneGSEA(TumorMacrophage, gene)

View(GSEAResults$sortkk)  

write.csv(GSEAResults$sortkk,"KEGG.csv")

GSEAResults2 <- GeneGSEA2(TumorMacrophage, gene)

View(GSEAResults2$sortkk) 

gseaplot2(GSEAResults$kk_gse, "hsa04510", color="firebrick", rel_heights=c(1,.2,.6), pvalue_table=T)
ggsave("kegg_1.pdf", width=8, height=6)

gseaplot2(GSEAResults$kk_gse, "hsa04151", color="firebrick", rel_heights=c(1,.2,.6), pvalue_table=T)
ggsave("kegg_2.pdf", width=8, height=6)

gseaplot2(GSEAResults$kk_gse, "hsa04390", color="firebrick", rel_heights=c(1,.2,.6), pvalue_table=T)
ggsave("kegg_3.pdf", width=8, height=6)

gseaplot2(GSEAResults$kk_gse, "hsa04520", color="firebrick", rel_heights=c(1,.2,.6), pvalue_table=T)
ggsave("kegg_4.pdf", width=8, height=6)

gseaplot2(GSEAResults$kk_gse, "hsa00010", color="firebrick", rel_heights=c(1,.2,.6), pvalue_table=T)
ggsave("kegg_5.pdf", width=8, height=6)

egmt <- GSEA(geneList, TERM2GENE=geneset, minGSSize=1, pvalueCutoff=0.9, verbose=FALSE)

gseaplot2(egmt, geneSetID=c(1), pvalue_table=T)

ggsave("kegg_Ferroptosis.pdf", width=8, height=6)






