#!/bin/sh
# sbatch --job-name=SSImp --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./sameSampleImputation_Cpp.sh
. /home/jrca253/anaconda2/etc/profile.d/conda.sh
conda activate r-environment
Rscript sameSampleImputation_Cpp.R
conda deactivate