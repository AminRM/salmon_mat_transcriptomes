#!/bin/bash
#SBATCH --job-name=hisat-mapping
#SBATCH --output dump.out
#SBATCH --error dump.err
#SBATCH --ntasks-per-node=10
#SBATCH --mem=20gb
#SBATCH --time=100:00:00
#SBATCH --mail-user=amin.esmail@csiro.au

module load bowtie
module load tophat
module load samtools

tophat2  --output-dir mapping /data/moh034/salmon/salmon_ref ${FILE}_R1.fastq.gz -2 ${FILE}_R2.fastq.gz

