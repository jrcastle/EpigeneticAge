#!/bin/sh
#sbatch --job-name=Train --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./TRAIN.sh

. /etc/profile.d/modules.sh
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "

# Modules needed for this R job
module load R/3.4.1

if [ -f glmnet.Training.RData ]; then rm glmnet.Training.RData; fi
if [ -f lambda.glmnet.Training.RData ]; then rm lambda.glmnet.Training.RData; fi
if [ -f MethAgevsSampleAge_trainsample.png ]; then rm MethAgevsSampleAge_trainsample.png; fi
if [ -f residual_hist_trainsample.png ]; then rm residual_hist_trainsample.png; fi

Rscript trainModel.R