#!/bin/sh
# sbatch --job-name=CONVERT --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./convertCHR-POStoCGNumber.sh

echo "python -u convertCHR-POStoCGNumber.py"
python -u convertCHR-POStoCGNumber.py

echo "Rscript cleanConvertedCHR-POStoCGNumber.R"
if [ -f ~/.Renviron ]; then mv ~/.Renviron ~/tmp.Renviron; fi
. /home/jrca253/anaconda2/etc/profile.d/conda.sh
conda activate r-environment
Rscript cleanConvertedCHR-POStoCGNumber.R
conda deactivate
if [ -f ~/tmp.Renviron ]; then mv ~/tmp.Renviron ~/.Renviron; fi
echo "rm tmp.csv"
rm tmp.csv
