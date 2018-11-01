#!/bin/sh
echo "python -u convertCHR-POStoCGNumber.py"
python -u convertCHR-POStoCGNumber.py
echo "Rscript cleanConvertedCHR-POStoCGNumber.R"
Rscript cleanConvertedCHR-POStoCGNumber.R
echo "rm tmp.csv"
rm tmp.csv