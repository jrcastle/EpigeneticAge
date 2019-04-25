#!/bin/sh
# sbatch --job-name=CLEAN --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./cleanConvertedCHR-POStoCGNumber.sh

if [ -f ~/.Renviron ]; then mv ~/.Renviron ~/tmp.Renviron; fi
. /home/jrca253/anaconda2/etc/profile.d/conda.sh
conda activate r-environment
Rscript cleanConvertedCHR-POStoCGNumber.R
conda deactivate
if [ -f ~/tmp.Renviron ]; then mv ~/tmp.Renviron ~/.Renviron; fi