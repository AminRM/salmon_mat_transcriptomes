#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output dump.out
#SBATCH --error dump.err
#SBATCH --ntasks-per-node=10
#SBATCH --mem=16gb
#SBATCH --time=10:00:00
#SBATCH --mail-user=amin.esmail@csiro.au

module load samtools/1.9.0

samtools merge ${FILE}_accepted_hits.bam ${FILE}_L001_output.bam ${FILE}_L002_output.bam
samtools sort -n -o ${FILE}_accepted_hits_sorted.bam ${FILE}_accepted_hits.bam
samtools view ${FILE}_accepted_hits_sorted.bam > ${FILE}_accepted_hits_sorted.sam
