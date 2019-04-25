#!/bin/sh
model_output="model_coefficients.csv"
out_file="ClockCpGs.txt"

if [ -f $out_file ]
then
    rm $out_file
fi

tr ',' '\t' < $model_output | awk {'print $1'} > tmp.txt
sed 's/\"//g' tmp.txt > tmp2.txt
sed '1d;2d' tmp2.txt > $out_file
rm tmp.txt
rm tmp2.txt