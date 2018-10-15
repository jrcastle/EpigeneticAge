#!/bin/sh
#sbatch --job-name=Train --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./MULTPLE_CV.sh

. /etc/profile.d/modules.sh
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "

# Modules needed for this R job
module load R/3.4.1

if [ -f df.lambda.RData ]; then rm df.lambda.RData; fi
if [ -f df.mse.RData ]; then rm df.mse.RData; fi
if [ -f lambda.min.RData ]; then rm lambda.min.RData; fi
if [ -f mse.min.RData ]; then rm mse.min.RData; fi

Rscript multipleCV.R