#!/bin/sh
# sbatch --job-name=IMPUTE --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./IMPUTE2.sh

. /etc/profile.d/modules.sh
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "

# Modules needed for this R job
module load R/3.4.1

#R Program execution command
echo 'ulimit -s unlimited'
ulimit -s unlimited
Rscript imputeMissingData2.R
