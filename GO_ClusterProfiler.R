if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
install.packages("tidyverse")
install.packages("Category")

# Load library and search for the organism
library(clusterProfiler)
library(pathview)
library(tidyverse)
library(enrichplot)
library(GenomicRanges)
library(Category)
library(AnnotationHub)
hub <- AnnotationHub()
yes
query(hub, c("salmo salar","orgdb"))
salmodb <- hub[["AH72154"]]
DatPkgFactory(salmodb)
columns(salmodb)

#read in expressed genes per tissue as background dataset 
ovary_gene_universe <- read.table("/Users/moh034/RNA_Seq/Ovary_gene_universe.txt", header=T, sep="\t")
pit_gene_universe <- read.table("/Users/moh034/RNA_Seq/pit_gene_universe.txt", header=T, sep="\t")
liver_gene_universe <- read.table("/Users/moh034/RNA_Seq/liver_gene_universe.txt", header=T, sep="\t")
head (ovary_gene_universe)
head (pit_gene_universe)
head (liver_gene_universe)

########################################################################################################
up- and downregulated gene clusters in ovary 
########################################################################################################
ovary_up <- read.table("/Users/moh034/RNA_Seq/cluster_1_ovary4clusterprofiler.txt", header=T, sep="\t")
dim(ovary_up)
head(ovary_up)
rownames(ovary_up)<- NULL
ovary_up_FC <- ovary_up$logFC
names(ovary_up_FC) <- ovary_up$geneID
head(ovary_up_FC)
ovary_up_FC <- sort(ovary_up_FC, decreasing = TRUE)
ego_ovary_up <- enrichGO(gene = names(ovary_up_FC), universe = as.character(ovary_gene_universe$ENTREZID), keyType = 'ENTREZID', OrgDb = salmodb, ont="ALL", pAdjustMethod = "fdr",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_ovary_up_Summary <- data.frame(ego_ovary_up)
write.csv(ego_ovary_up_Summary, file = "/Users/moh034/RNA_Seq/ego_ovary_up.csv", row.names=FALSE)
dotplot(ego_ovary_up, showCategory=20, font.size = 10)
dev.off()

#ovary T4 downregulated  cluster
ovary_down <- read.table("/Users/moh034/RNA_Seq/cluster_2_ovary4clusterprofiler.txt", header=T, sep="\t")
dim(ovary_down)
head(ovary_down)
rownames(ovary_downT4)<- NULL
ovary_down_FC <- ovary_down$logFC
names(ovary_down_FC) <- ovary_down$geneID
head(ovary_down_FC)
ovary_down_FC <- sort(ovary_down_FC, decreasing = TRUE)
ego_ovary_down <- enrichGO(gene = names(ovary_down_FC), universe = as.character(ovary_gene_universe$ENTREZID), keyType = 'ENTREZID', OrgDb = salmodb, ont="ALL", pAdjustMethod = "fdr",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_ovary_down_Summary <- data.frame(ego_ovary_down)
write.csv(ego_ovary_down_Summary, file = "/Users/moh034/RNA_Seq/ego_ovary_down.csv", row.names=FALSE)
dotplot(ego_ovary_down, showCategory=20, font.size = 10)
dev.off()

########################################################################################################
up- and downregulated clusters in pituitary  
########################################################################################################
#pituitary upregulated  cluster
pit_up <- read.table("/Users/moh034/RNA_Seq/Pit_Up_Clust4ClusterProfiler.txt", header=T, sep="\t")
dim(pit_up)
head(pit_up)
rownames(pit_up)<- NULL
pit_up_FC <- pit_up$logFC
names(pit_up_FC) <- pit_up$geneID
head(pit_up_FC)
pit_up_FC <- sort(pit_up_FC, decreasing = TRUE)
ego_pit_up <- enrichGO(gene = names(pit_up_FC), universe = as.character(pit_gene_universe$ENTREZID), keyType = 'ENTREZID', OrgDb = salmodb, ont="ALL", pAdjustMethod = "fdr",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_pit_up_Summary <- data.frame(ego_pit_up)
write.csv(ego_pit_up_Summary, file = "/Users/moh034/RNA_Seq/ego_pit_up.csv", row.names=FALSE)
dotplot(ego_pit_up, showCategory=20, font.size = 10)
dev.off()

#pituitary downregulated  cluster
pit_down <- read.table("/Users/moh034/RNA_Seq/Pit_DownClust4clusterProfiler.txt", header=T, sep="\t")
dim(pit_down)
head(pit_down)
rownames(pit_down)<- NULL
pit_down_FC <- pit_down$logFC
names(pit_down_FC) <- pit_down$geneID
head(pit_down_FC)
pit_down_FC <- sort(pit_down_FC, decreasing = TRUE)
ego_pit_down <- enrichGO(gene = names(pit_down_FC), universe = as.character(pit_gene_universe$ENTREZID), keyType = 'ENTREZID', OrgDb = salmodb, ont="ALL", pAdjustMethod = "fdr",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_pit_down_Summary <- data.frame(ego_pit_down)
write.csv(ego_pit_down_Summary, file = "/Users/moh034/RNA_Seq/ego_pit_down.csv", row.names=FALSE)
dotplot(ego_pit_down, showCategory=20, font.size = 10)
dev.off()
