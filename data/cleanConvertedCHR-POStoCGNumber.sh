#!/bin/sh
# sbatch --job-name=CLEAN --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./cleanConvertedCHR-POStoCGNumber.sh

. /home/jrca253/anaconda2/etc/profile.d/conda.sh
conda activate r-environment
Rscript cleanConvertedCHR-POStoCGNumber.R
conda deactivate