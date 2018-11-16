#!/bin/sh
# sbatch --job-name=KeepCpGs --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./KEEP_SPECIFIC_CPGS.sh
. /etc/profile.d/modules.sh
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "

# Modules needed for this R job
module load R/3.4.1

#R Program execution command
Rscript keep_specific_cpgs.R