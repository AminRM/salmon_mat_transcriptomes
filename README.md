# salmon_maturation_transcriptomics
Mohamed et al 2020 (under review) Integrated transcriptome, DNA methylome and chromatin state accessibility landscapes reveal regulators of Atlantic salmon maturation

This repository contains scripts to perform RNA-seq analysis from 3 salmon tissues (pituitary gland, ovary and liver) aiming at characetrizing transcriptional remodelling associated with onset of puberty in salmon. This will identify differentially expressed genes and clusters of co-expressed genes.  

the following scripts were used to perform dataQC, map RNA-seq data to the Atlantic salmon reference genome ICSASG_v2 (Lien et al., 2016), perform explorartory analyese, perform differential expression (DE) analyses using the R package edgeR, perfom co-expression analyses

Raw data QC: fastQC.sh 

genome mapping: bowtie2-build.sh and tophat2.sh were used to index the salmon genome and to map ATAC-seq data to the salmon genome using tophat2 

Read counting using HTSeq_count.sh to get raw counts before merge_sort_bam2sam.sh was used to merge/sort bam files per sample and convert BAM to SAM format

build gene expression "raw counts" matrix using HTSeq2CountMatrix.R and perform data QC to check samples agreement, detect outliers and perform DE analyses  using Explor_DE_Analysis 

functional profile : GO_ClusterProfiler.R will perform GO enrichment analyese on the list oc co-expressed clusters using ClusterProfiler 
