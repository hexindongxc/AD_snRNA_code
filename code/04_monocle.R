
library(monocle)
table(scRNA_celltype$minor)
Macro <- scRNA_celltype
sample_ann <-  Macro@meta.data 

gene_ann <- data.frame(gene_short_name = rownames(Macro@assays$RNA),
                       row.names =  rownames(Macro@assays$RNA))
pd <- new("AnnotatedDataFrame",data=sample_ann)
fd <- new("AnnotatedDataFrame",data=gene_ann)
ct=as.data.frame(Macro@assays$RNA$counts)
sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F)
cds <- reduceDimension(cds, max_components = 2, num_dim = 25,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 5) 
plot_cell_clusters(cds, 1, 2 )
colnames(pData(cds))
pData(cds)$Cluster=pData(cds)$minor
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~Cluster")
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Cluster")
ggsave("clusters.pdf",width=10,height=10)
plot_cell_trajectory(cds, color_by = "seurat_clusters")
ggsave("clusters.pdf",width=10,height=10)
plot_cell_trajectory(cds, color_by = "Pseudotime") + theme(legend.position = "right")
ggsave("Pseudotime.pdf",width=8,height=6)
plot_cell_trajectory(cds, color_by = "group")+ theme(legend.position = "right")
ggsave("type.pdf",width=8,height=6)

plot_cell_trajectory(cds, color_by = "State")

ggsave("State.pdf",width=10,height=10)

plot_cell_trajectory(cds, color_by = "minor")+ theme(legend.position = "right")
ggsave("Minor.pdf",width=8,height=6)

blast_genes <- row.names(subset(fData(cds),
                                gene_short_name %in% c("RRBP1")))
plot_genes_jitter(cds[blast_genes,],
                  grouping = "Pseudotime",
                  min_expr = 0.1)

saveRDS(cds, "monocle2.rds")
table(Macro$seurat_clusters)
cds <- Macro
p1 <-plot_cell_trajectory(cds, color_by = "minor")
p2 <- plot_cell_trajectory(cds, color_by = "seurat_clusters")
p1+p2
p1
ggsave("macro_type.pdf",p1,width=10,height=10)

library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "macro")  + scale_color_npg() 
p2=plot_cell_trajectory(cds, color_by = "State")  + scale_color_nejm()
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
p3=plot_cell_trajectory(cds, color_by = "State")  + scale_color_manual(values = colour)
p1|p2|p3

p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "macro") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = col_vector) 
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "macro")+
  scale_color_manual(values = col_vector) +
  theme(legend.title = element_blank()) 
p1|p2
ggsave("macro_monocle.pdf",width=16,height=8)

p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "group") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = col_vector) 
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "group")+
  scale_color_manual(values = col_vector) +
  theme(legend.title = element_blank()) 
p1|p2
ggsave("group_monocle.pdf",width=16,height=8)

p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "seurat_clusters") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = col_vector) 
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "seurat_clusters")+
  scale_color_manual(values = col_vector) +
  theme(legend.title = element_blank()) 
p1|p2
ggsave("cluster_monocle.pdf",width=16,height=8)

library(ggplot2)
df <- pData(cds) 
ggplot(df, aes(Pseudotime, colour = minor, fill=minor)) +
  geom_density(bw=0.5,linewidth=1,alpha = 0.5)
ggsave("minor_monocle.pdf",width=10,height=8)

ggplot(df, aes(Pseudotime, colour = group, fill=group)) +
  geom_density(bw=0.5,linewidth=1,alpha = 0.5)
ggsave("minor_monocle.pdf",width=10,height=8)


diff.genes <- FindAllMarkers(scRNA_celltype, logfc.threshold=0.5,
               test.use="wilcox", 
               return.thresh=0.25, 
               min.pct=0.3, only.pos=F)

sig_diff.genes <- subset(diff.genes,p_val_adj<0.00001&abs(avg_log2FC)>0.85)$gene
sig_diff.genes <- unique(as.character(sig_diff.genes))

diff_test <- differentialGeneTest(cds[sig_diff.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))

write.table(sig_gene_names, file="gene_monocle.txt",quote = F,
            sep = "\t",row.names = F)
p1 = plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters=3,
                             show_rownames=F, return_heatmap=T)
p1
ggsave("pseudotime_heatmap1.pdf", plot = p1, width = 10, height = 10)

disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 2&dispersion_empirical >= 0.9*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(cds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-02))
sig_gene_names <- row.names(diff_test)
p2 = plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters=3,
                             show_rownames=T, return_heatmap=T)
p2
ggsave("pseudotime_heatmap.pdf", plot = p2, width = 6, height = 16)
write.table(sig_gene_names, file="gene_monocle.txt",quote = F,
            sep = "\t",row.names = F)

  s.genes <- c("RRBP1","PLEC","MUC5B","PLCG2","RNF213","PPP1R15A")
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "group", color_by = "group")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "group", color_by = "group")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "minor")
plotc <- p1|p2|p3
plotc
ggsave("Genes_Jitterplot.pdf",  width = 10, height = 16)

plot_genes_in_pseudotime(cds[s.genes,], color_by = "group")


