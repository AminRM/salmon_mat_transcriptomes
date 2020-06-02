#We have 64 samples (4 tissue "B, P, O, L" x 4 time points "T1, T2, T3, T4" x 4 reps "F1, F2, F3, F4")
#sample names "header of the count matrix" will be needed at the end 

samples <- c("B_T1F1","B_T1F2","B_T1F3","B_T1F4","B_T2F1","B_T2F2","B_T2F3","B_T2F4","B_T3F1","B_T3F2","B_T3F3","B_T3F4","B_T4F1","B_T4F2","B_T4F3","B_T4F4", 
             "P_T1F1","P_T1F2","P_T1F3","P_T1F4","P_T2F1","P_T2F2","P_T2F3","P_T2F4","P_T3F1","P_T3F2","P_T3F3","P_T3F4","P_T4F1","P_T4F2","P_T4F3","P_T4F4",
             "O_T1F1","O_T1F2","O_T1F3","O_T1F4","O_T2F1","O_T2F2","O_T2F3","O_T2F4","O_T3F1","O_T3F2","O_T3F3","O_T3F4","O_T4F1","O_T4F2","O_T4F3","O_T4F4",
             "L_T1F1","L_T1F2","L_T1F3","L_T1F4","L_T2F1","L_T2F2","L_T2F3","L_T2F4","L_T3F1","L_T3F2","L_T3F3","L_T3F4","L_T4F1","L_T4F2","L_T4F3","L_T4F4")

#read samples function 
read.sample <- function(sample.name) {
  file.name <- paste(sample.name, ".counts", sep="")
  result <- read.delim(file.name, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"))
}
sample.1 <- read.sample(samples[1])
sample.2 <- read.sample(samples[2])
sample.3 <- read.sample(samples[3])
sample.4 <- read.sample(samples[4])
sample.5 <- read.sample(samples[5])
sample.6 <- read.sample(samples[6])
sample.7 <- read.sample(samples[7])
sample.8 <- read.sample(samples[8])
sample.9 <- read.sample(samples[9])
sample.10 <- read.sample(samples[10])
sample.11 <- read.sample(samples[11])
sample.12 <- read.sample(samples[12])
sample.13 <- read.sample(samples[13])
sample.14 <- read.sample(samples[14])
sample.15 <- read.sample(samples[15])
sample.16 <- read.sample(samples[16])
sample.17 <- read.sample(samples[17])
sample.18 <- read.sample(samples[18])
sample.19 <- read.sample(samples[19])
sample.20 <- read.sample(samples[20])
sample.21 <- read.sample(samples[21])
sample.22 <- read.sample(samples[22])
sample.23 <- read.sample(samples[23])
sample.24 <- read.sample(samples[24])
sample.25 <- read.sample(samples[25])
sample.26 <- read.sample(samples[26])
sample.27 <- read.sample(samples[27])
sample.28 <- read.sample(samples[28])
sample.29 <- read.sample(samples[29])
sample.30 <- read.sample(samples[30])
sample.31 <- read.sample(samples[31])
sample.32 <- read.sample(samples[32])
sample.33 <- read.sample(samples[33])
sample.34 <- read.sample(samples[34])
sample.35 <- read.sample(samples[35])
sample.36 <- read.sample(samples[36])
sample.37 <- read.sample(samples[37])
sample.38 <- read.sample(samples[38])
sample.39 <- read.sample(samples[39])
sample.40 <- read.sample(samples[40])
sample.41 <- read.sample(samples[41])
sample.42 <- read.sample(samples[42])
sample.43 <- read.sample(samples[43])
sample.44 <- read.sample(samples[44])
sample.45 <- read.sample(samples[45])
sample.46 <- read.sample(samples[46])
sample.47 <- read.sample(samples[47])
sample.48 <- read.sample(samples[48])
sample.49 <- read.sample(samples[49])
sample.50 <- read.sample(samples[50])
sample.51 <- read.sample(samples[51])
sample.52 <- read.sample(samples[52])
sample.53 <- read.sample(samples[53])
sample.54 <- read.sample(samples[54])
sample.55 <- read.sample(samples[55])
sample.56 <- read.sample(samples[56])
sample.57 <- read.sample(samples[57])
sample.58 <- read.sample(samples[58])
sample.59 <- read.sample(samples[59])
sample.60 <- read.sample(samples[60])
sample.61 <- read.sample(samples[61])
sample.62 <- read.sample(samples[62])
sample.63 <- read.sample(samples[63])
sample.64 <- read.sample(samples[64])

nrow(sample.1) == nrow(sample.2)
all(sample.1$gene == sample.2$gene)

#add geneIDs as first column and then add the counts for each sample
64RawCount.data <- sample.1
64RawCount.data <- cbind(sample.1, sample.2$count, sample.3$count, sample.4$count, sample.5$count, sample.6$count, sample.7$count,sample.8$count,sample.9$count, sample.10$count, sample.11$count, sample.12$count, sample.13$count, sample.14$count,sample.15$count, sample.16$count,
                  sample.17, sample.18$count, sample.19$count, sample.20$count, sample.21$count, sample.22$count, sample.23$count,sample.24$count,sample.25$count, sample.26$count, sample.27$count, sample.28$count, sample.29$count, sample.30$count,sample.31$count, sample.32$count,
                  sample.33, sample.34$count, sample.35$count, sample.36$count, sample.37$count, sample.38$count, sample.39$count,sample.40$count,sample.41$count, sample.42$count, sample.43$count, sample.44$count, sample.45$count, sample.46$count,sample.47$count, sample.48$count,
                  sample.49, sample.50$count, sample.51$count, sample.52$count, sample.53$count, sample.54$count, sample.55$count,sample.56$count,sample.57$count, sample.58$count, sample.59$count, sample.60$count, sample.61$count, sample.62$count,sample.63$count, sample.64$count
                 )
#take away last 5 lines from HTSeq count files 
RawCount.data <- all.data[1:(nrow(all.data)-5),]
#add column names 
colnames(64RawCount.data) <- samples
#Final check!
dim(RawCount.data)
head(RawCount.data)
tail(RawCount.data)

#read in a file that contains gene length; this will be used later to get FPKM values 
file_len<- read.delim("gene_length.txt",header=F,sep="\t")


