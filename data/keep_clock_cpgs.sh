#!/bin/sh
cpg_file="tmp.txt"
meth_file="meth_N_lt10_missing.txt"
out_file="tmp2.txt"

# IF OUTFILE EXISTS, REMOVE
if [ -f $out_file ]
then
    rm $out_file
fi

# HEADER
head -1 $meth_file > $out_file

# CpGs
i=1
for cpg in `cat $cpg_file`; do
    echo "Finding CpG $i ..."
    grep $cpg $meth_file >> $out_file
    let "i++"
done

echo 'DONE!'

