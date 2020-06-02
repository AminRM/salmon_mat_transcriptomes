#! /bin/bash
#SBATCH --output  dump.out
#SBATCH --error dump.err
#SBATCH --mem=5gb
#SBATCH --time=30:00
#SBATCH --mail-user=amin.esmail@csiro.au
/apps/fastqc/0.11.8/fastqc ${FILE}
