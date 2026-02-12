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
seurat_data <- read.csv("./GSE138852_counts.csv", row.names = 1, header = TRUE)
seurat_obj <- CreateSeuratObject(counts = seurat_data, min.features = 200, min.cells = 3)
sce <- seurat_obj
sce@meta.data$orig.ident <- gsub("^[^_]*_", "", rownames(sce@meta.data))
table(sce$orig.ident)
sce@meta.data$group <- NA
sce@meta.data$group <- ifelse(sce@meta.data$orig.ident%in%c("AD1_AD2","AD3_AD4","AD5_AD6"),"AD",sce@meta.data$group)
sce@meta.data$group <- ifelse(sce@meta.data$orig.ident%in%c("Ct1_Ct2","Ct3_Ct4","Ct5_Ct6"),"CT",sce@meta.data$group)
table(sce@meta.data$group)
#AD   CT 
#6673 6541 
sce[["percent.ercc"]] <- PercentageFeatureSet(sce, pattern="^ERCC-")  
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern="^MT-")   
sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern="^RP[SL]")  
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(sce))
sce[["percent.HB"]] <- PercentageFeatureSet(sce, features=HB.genes)  
raw_cell <- sce@meta.data
raw_count <- table(raw_cell$orig.ident)
raw_count
sce$orig.ident <- factor(sce$orig.ident,levels=c("Ct1_Ct2","Ct3_Ct4","Ct5_Ct6","AD1_AD2","AD3_AD4","AD5_AD6"))
pearplot_befor <- VlnPlot(sce, group.by="orig.ident",
                          features=c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                          pt.size=0, ncol=4)
pearplot_befor
ggsave(pearplot_befor, file="pearplot_befor.pdf",  width=12, height=6)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category=="qual",]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
sample_color <- col_vector[1:33]
sample_color
Feature_ber1 <- FeatureScatter(sce, feature1="nFeature_RNA",
                               feature2="nCount_RNA",
                               group.by="orig.ident",
                               cols=sample_color)
Feature_ber2 <- FeatureScatter(sce, feature1="percent.mt",
                               feature2="nCount_RNA",
                               group.by="orig.ident",
                               cols=sample_color)
Feature_ber3 <- FeatureScatter(sce, feature1="percent.mt",
                               feature2="nFeature_RNA",
                               group.by="orig.ident",
                               cols=sample_color)
Feature_ber1 <- Feature_ber1 + theme(legend.position="none")
Feature_ber2 <- Feature_ber2 + theme(legend.position="none")
Feature_ber <- ggarrange(Feature_ber1, Feature_ber2, Feature_ber3, ncol=3, nrow=1, widths=c(1,1,1.2))
Feature_ber
ggsave(Feature_ber, file="Feature_before.pdf", width=12, height=6)
sce <- subset(sce, subset = nFeature_RNA<1600 & nCount_RNA<3000 & percent.mt<5 &percent.Ribo<5)
pearplot_after <- VlnPlot(sce, group.by="orig.ident",
                          features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.Ribo"),
                          pt.size=0,
                          ncol=4)
pearplot_after
ggsave(pearplot_after, file="pearplot_after.pdf", width=12, height=6)
Feature_ber1 <- FeatureScatter(sce, feature1="nFeature_RNA",
                               feature2="nCount_RNA",
                               group.by="orig.ident",
                               cols=sample_color)
Feature_ber2 <- FeatureScatter(sce, feature1="percent.mt",
                               feature2="nCount_RNA",
                               group.by="orig.ident",
                               cols=sample_color)
Feature_ber3 <- FeatureScatter(sce, feature1="percent.mt",
                               feature2="nFeature_RNA",
                               group.by="orig.ident",
                               cols=sample_color)
Feature_ber1 <- Feature_ber1 + theme(legend.position="none")
Feature_ber2 <- Feature_ber2 + theme(legend.position="none")
Feature_ber <- ggarrange(Feature_ber1, Feature_ber2, Feature_ber3, ncol=3, nrow=1, widths=c(1,1,1.2))
Feature_ber
ggsave(Feature_ber, file="Feature_after.pdf", width=12, height=6)

clean_cell <- sce@meta.data
clean_count <- table(clean_cell$orig.ident)
clean_count
sum(clean_count) 

sce_filtered <- JoinLayers(sce)  
sce_processed <- NormalizeData(sce_filtered) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

library(DoubletFinder)

if (!exists("paramSweep_v3")) {
  paramSweep_function <- paramSweep
  doubletFinder_function <- doubletFinder
} else {
  paramSweep_function <- paramSweep_v3
  doubletFinder_function <- doubletFinder_v3
}

find_doublets <- function(seurat_obj, sample_name) {
  subset_obj <- subset(seurat_obj, subset = orig.ident == sample_name)
  subset_obj <- NormalizeData(subset_obj) %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(verbose = FALSE)
  sweep.res <- paramSweep_function(subset_obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal.pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cat("样本", sample_name, "的最佳pK值为:", optimal.pk, "\n")
  nExp_poi <- round(0.075 * ncol(subset_obj))
  subset_obj <- doubletFinder_function(subset_obj, 
                                       PCs = 1:20, 
                                       pN = 0.25, 
                                       pK = optimal.pk, 
                                       nExp = nExp_poi, 
                                       reuse.pANN = FALSE, 
                                       sct = FALSE)
  
  return(subset_obj)
}

doublet_results <- list()
for (sample in unique(sce_processed$orig.ident)) {
  doublet_results[[sample]] <- find_doublets(sce_processed, sample)
}

doublet_class <- c()
for (sample in names(doublet_results)) {
  sample_obj <- doublet_results[[sample]]
  df_column <- grep("DF.classifications", names(sample_obj@meta.data), value = TRUE)
  sample_cells <- colnames(sample_obj)
  sample_class <- sample_obj@meta.data[[df_column]]
  names(sample_class) <- sample_cells
  doublet_class <- c(doublet_class, sample_class)
}

sce_processed$doublet_class <- doublet_class[colnames(sce_processed)]
doublet_summary <- table(sce_processed$orig.ident, sce_processed$doublet_class)
print(doublet_summary)

p10 <- DimPlot(sce_processed, reduction = "pca", group.by = "doublet_class", split.by = "orig.ident")
print(p10)
ggsave("05_doublet_detection.pdf", p10, width = 12, height = 5)

sce_clean <- subset(sce_processed, subset = doublet_class == "Singlet")


scRNA_harmony <- NormalizeData(sce_clean) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose=FALSE)
plot1 <- VariableFeaturePlot(scRNA_harmony)
plot2 <- LabelPoints(plot=plot1, points=head(VariableFeatures(scRNA_harmony),10), repel=TRUE, size=2.5)
plot2
ggsave(plot2, file="VariableFeatures.pdf", width=6, height=6)
DimPlot(scRNA_harmony, reduction="pca", group.by="orig.ident")
ggsave(file="PCA_Dimplot.pdf", width=8, height=6)
elbowplot <- ElbowPlot(scRNA_harmony, ndims=50)  
elbowplot
ggsave(elbowplot, file="elbowplot.pdf", width=6, height=6)
pc.nums=25
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
clustreeplot <- clustree(scRNA_harmony@meta.data, prefix="RNA_snn_res.")  #0.2
clustreeplot
ggsave(clustreeplot, file="clustreeplot.pdf", width=12, height=10)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction="harmony", dims=1:pc.nums)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution=0.2) ##
table(scRNA_harmony@active.ident)  
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction="harmony", dims=1:pc.nums)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction="harmony", dims=1:pc.nums)
plot1 <- DimPlot(scRNA_harmony, reduction="umap", group.by="orig.ident")
plot2 <- DimPlot(scRNA_harmony, reduction="umap", label=T)
plot1
plot2
ggsave(plot1, file="umapplot1.pdf", width=8, height=6)
ggsave(plot2, file="umapplot2.pdf", width=8, height=6)
saveRDS(scRNA_harmony, "scRNA_harmony.rds")
