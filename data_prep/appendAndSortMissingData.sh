#!/bin/sh
# sbatch --job-name=APP_SORT --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./appendAndSortMissingData.sh

if [ -f ~/.Renviron ]; then mv ~/.Renviron ~/tmp.Renviron; fi
. /home/jrca253/anaconda2/etc/profile.d/conda.sh
conda activate r-environment
Rscript appendAndSortMissingData.R
conda deactivate
if [ -f ~/tmp.Renviron ]; then mv ~/tmp.Renviron ~/.Renviron; fi
