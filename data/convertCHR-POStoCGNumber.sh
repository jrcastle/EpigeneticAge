#!/bin/sh
# sbatch --job-name=CONVERT --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./convertCHR-POStoCGNumber.sh

echo "python -u convertCHR-POStoCGNumber.py"
python -u convertCHR-POStoCGNumber.py

echo "Rscript cleanConvertedCHR-POStoCGNumber.R"
. /home/jrca253/anaconda2/etc/profile.d/conda.sh
conda activate r-environment
Rscript cleanConvertedCHR-POStoCGNumber.R
conda deactivate

echo "rm tmp.csv"
rm tmp.csv