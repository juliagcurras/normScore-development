#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 00:01:00
#SBATCH --mem-per-cpu  1GB
#SBATCH --error=/home/ulc/es/jgc/Simulacions/error_item0156.txt
#SBATCH --output=/mnt/lustre/scratch/nlsas/home/ulc/es/jgc/Simulations/output_item0156.log

module load cesga/system R/4.4.2

export R_LIBS_USER=$STORE/Rlibs/4.4.2

Rscript execute2.R