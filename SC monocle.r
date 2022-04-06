
library(Seurat);library(ggsci);library(ggplot2);library(cowplot);library(dplyr);library(clustree)
SC <- readRDS('/home/data/small_wound/SC/SCs in small_wound_Seurat.rds')
library(monocle)
wa<- SC
data <- as(as.matrix(wa@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = wa@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

HSMM<-monocle_cds

##
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
print(head(pData(HSMM)))
######################################
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
pdf("choosing genes.pdf")
plot_ordering_genes(HSMM)
dev.off()
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)
pdf("cell_trajectory state.pdf")
plot_cell_trajectory(HSMM, color_by = "State",cell_size=)
dev.off()
pdf("cell_trajectory pseudotime.pdf")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()
head(pData(HSMM))
pdf("cell_trajectory by time.pdf")
plot_cell_trajectory(HSMM, color_by = "wall")
dev.off() 


########### Finding Genes that Change as a Function of Pseudotime ###############
my_cds_subset<- HSMM
my_pseudotime_de <- differentialGeneTest(my_cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 4)
my_pseudotime_de %>% arrange(qval) %>% head()%>% select(gene_short_name) -> my_pseudotime_gene
write.csv(my_pseudotime_gene, "my_pseudotime_gene.csv")
my_pseudotime_gene <- my_pseudotime_gene$gene_short_name
my_pseudotime_gene

pdf(file="10_Pseudotimegenes_in_pseudotime.pdf")
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])
dev.off()

##### Clustering Genes by Pseudotemporal Expression Pattern ##############

my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
pdf("13_pseudotime_heatmap.pdf")
my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,], num_clusters = 3, cores = 1, show_rownames = TRUE, return_heatmap = TRUE)
dev.off()

newdata <- data.frame(Pseudotime = seq(min(pData(my_cds_subset)$Pseudotime), max(pData(my_cds_subset)$Pseudotime), length.out = 100))
my_cluster <- cutree(my_pseudotime_cluster$tree_row, 3)
write.csv(my_cluster, file="13_pseudotime_cluster.csv")


########## Analyzing Branches in Single-Cell Trajectories ####################

BEAM_res <- BEAM(my_cds_subset, cores = 6)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
pdf("14_branch_heatmap.pdf")
my_branched_heatmap<- plot_genes_branched_heatmap(my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),], branch_point = 1, num_clusters = 4, cores = 1, use_gene_short_name = TRUE, show_rownames = TRUE, return_heatmap = TRUE)
dev.off()
dim(BEAM_res)

pdf("14_branch_heatmap_top50.pdf")
my_branched_heatmap<- plot_genes_branched_heatmap(my_cds_subset[row.names(BEAM_res[1:50,]),], branch_point = 1, num_clusters = 2, cores = 2, use_gene_short_name = TRUE, show_rownames = TRUE, return_heatmap = TRUE)
dev.off()
write.csv(BEAM_res, file="15_branch_gene.csv")

my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster, gene = row.names(my_row), stringsAsFactors = FALSE)
my_gene <- row.names(subset(fData(my_cds_subset), gene_short_name %in% head(my_row[my_row$cluster == 3,'gene'])))
#################################################################################
######################################################################
SC.markers <- FindAllMarkers( SC ,only.pos = TRUE, min.pct = 0.25, log2fc.threshold =0.25)
top10 <- SC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

library(celltalker)
set.seed(02221989)
head(ramilowski_pairs)
ligs <- as.character(unique(ramilowski_pairs$ligand))
ligs
 
library(biomaRt)
human=useMart('ensembl',dataset='hsapiens_gene_ensembl')
mouse=useMart('ensembl',dataset='mmusculus_gene_ensembl')
genes =ligs
g=getLDS(attributes=c("hgnc_symbol"),filters="hgnc_symbol",values=genes,mart=human,
     attributesL=c("mgi_symbol"),martL=mouse,uniqueRows=T)

ligs_list <- as.character(g[,2])



sigene <-rownames(my_pseudotime_gene)
sig_lig <-sigene[sigene%in%ligs_list]
sig_genes[sig_lig,]%>% arrange(qval) %>% head(50) -> sig_lig_res
pdf("16_significant_ligs_in_pseudotime_heatmap.pdf")
plot_pseudotime_heatmap(HSMM[sig_lig,],show_rownames=TRUE)
dev.off()
















