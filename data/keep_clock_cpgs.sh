#!/bin/sh
cpg_file="ClockCpGs.txt"
meth_file_K="meth_K_cpgs_in_KNT_imputed_vali.txt"
out_file_K="meth_K_cpgs_in_KNT_imputed_vali_ClockCpGs.txt"
meth_file_N="meth_N_cpgs_in_KNT_imputed.txt"
out_file_N="meth_N_cpgs_in_KNT_imputed_ClockCpGs.txt"
meth_file_T="meth_T_cpgs_in_KNT_imputed.txt"
out_file_T="meth_T_cpgs_in_KNT_imputed_ClockCpGs.txt"


# IF OUTFILE EXISTS, REMOVE
if [ -f $out_file_K ]
then
    rm $out_file_K
fi

if [ -f $out_file_N ]
then
    rm $out_file_N
fi

if [ -f $out_file_T ]
then
    rm $out_file_T
fi


# HEADER
head -1 $meth_file_K > $out_file_K
head -1 $meth_file_N > $out_file_N
head -1 $meth_file_T > $out_file_T


# CpGs
i=1
for cpg in `cat $cpg_file`; do
    echo "Finding CpG $i ..."
    grep -w $cpg $meth_file_K >> $out_file_K
    grep -w $cpg $meth_file_N >> $out_file_N
    grep -w $cpg $meth_file_T >> $out_file_T
    let "i++"
done

echo 'DONE!'

