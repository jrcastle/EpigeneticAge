#!/bin/sh
# sbatch --job-name=SSImp --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./sameSampleImputation_Cpp.sh

if [ -f ~/.Renviron ]; then mv ~/.Renviron ~/tmp.Renviron; fi
. /home/jrca253/anaconda2/etc/profile.d/conda.sh
conda activate r-environment
Rscript sameSampleImputation_Cpp.R
conda deactivate
if [ -f ~/tmp.Renviron ]; then mv ~/tmp.Renviron ~/.Renviron; fi
