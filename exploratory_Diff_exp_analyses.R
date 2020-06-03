#starting from the matrix containing the raw count data 
#using utilities "functions" within trinity to explore your data

1- perfom log2 and CPM transformation to obtain multi-pages results of number of mapped fragments per treatment and samples correlation 
providing a samples.txt file describing treatment-replicate relationships 

module load R 
/apps/trinity/2.8.4/Analysis/DifferentialExpression/PtR --matrix RawCount.data --samples samples.txt --compare_replicates --CPM --log2 --min_rowSums 10

2-Obtain TMM normalised-FPKM matrix using edgeR 
library(edgeR)
rnaseqMatrix = read.table("RawCount.data", header=T, row.names=1, com='', check.names=F)
rnaseqMatrix = as.matrix(rnaseqMatrix)
rnaseqMatrix = round(rnaseqMatrix)
exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study = calcNormFactors(exp_study) #the default method for computing scaling factors uses a trimmed mean of M-values (TMM) between each pair of samples. 
exp_study$samples$eff.lib.size = exp_study$samples$lib.size * exp_study$samples$norm.factors
logCPM <- cpm(exp_study, prior.count=2, log=TRUE)
FPKM <- logCPM(y) #where y is the gene length infor  
write.table(FPKM, file="FPKMCount.data", quote=F, sep="\t", row.names=F)

3- check your samples for outliers, compare replicates across all samples by generating a correlation matrix for all sample replicates using log2-FPKM (TMM normalised matrix)
module load R 
/apps/trinity/2.8.4/Analysis/DifferentialExpression/PtR --matrix FPKMCount.data --min_rowSums 10 -s samples.txt --sample_cor_matrix
#samples separated by tissues and T1 samples (control) clearly separted from post-maturation(T2-T4) samples 

4- DE analysis liver T2 vs T1, I will annotate the first DE test only as the rest will be very similar: 
this analyses has been repeated for all comparisons T2, T3, T4 vs T1 at each tissue, the matrix has all expression data for 64 samples 

#Read the raw counts matrix and choose samples and perform 
data = read.table("/flush1/moh034/Maturation/all_data_matrix/all_data_new_matrix.matrix", header=T, row.names=1, com='')
col_ordering = c(5,6,7,8,1,2,3,4)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
#filter lowly expressed genes
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
#define the replicates for each condition 
conditions = factor(c(rep("Liver_T2", 4), rep("Liver_T1", 4)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
#test for DE 
et = exactTest(exp_study, pair=c("Liver_T2", "Liver_T1"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="Liver_T2", sampleB="Liver_T1", result_table)
result_table$logFC = -1 * result_table$logFC
#write the results table into a txt file 
write.table(result_table, file='Liver_T2_vs_Liver_T1.DE_results', sep='	', quote=F, row.names=T)
#write the matrix of the raw counts of these genes into a txt file
write.table(rnaseqMatrix, file='Liver_T2_vs_Liver_T1.count_matrix', sep='	', quote=F, row.names=T)
#plot results as MA plot!
pdf("Liver_T2_vs_Liver_T1.edgeR.DE_results.MA_n_Volcano.pdf")
source("/flush1/moh034/Maturation/rnaseq_plot_funcs.R")
plot_MA(rownames(result_table), result_table$logCPM, result_table$logFC, result_table$FDR)
dev.off()







