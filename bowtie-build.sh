#!/bin/bash
#SBATCH --job-name=tophat-index
#SBATCH --output dump.out
#SBATCH --error dump.err
#SBATCH --mem=20gb
#SBATCH --time=10:00:00
#SBATCH --mail-user=amin.esmail@csiro.au

module load bowtie
module load tophat

bowtie2-build GCA_000233375.4_ICSASG_v2_genomic.fna salmon_ref








