#!/bin/bash
#SBATCH --job-name=br_counts
#SBATCH --time=20:00:00
#SBATCH --output dump.out
#SBATCH --error dump.err
#SBATCH --mem=10gb
#SBATCH --mail-user=amin.esmail@csiro.au
module load python
python -m HTSeq.scripts.count -a 10 --stranded=reverse -r name ${FILE}_accepted_hits_sorted.sam /data/moh034/salmon/Salmo_salar-annotation.gff3 > ${FILE}.counts

