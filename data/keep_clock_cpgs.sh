#!/bin/sh
cpg_file="ClockCpGs.txt"
meth_file="meth_K_cpgs_in_KNT_imputed_vali.txt"
out_file="meth_K_cpgs_in_KNT_imputed_vali_ClockCpGs.txt"

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
    grep -w $cpg $meth_file >> $out_file
    let "i++"
done

echo 'DONE!'

