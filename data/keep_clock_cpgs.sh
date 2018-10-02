#!/bin/sh
OUT="meth_CGNumber_ClockCpGs.txt"

# IF OUTFILE EXISTS, REMOVE
if [ -f $OUT ]
then
    rm $OUT
fi

# HEADER
head -1 meth_CGNumber.txt > $OUT

# CpGs
i=1
for cpg in `cat ClockCpGs.txt`; do
    echo "Finding CpG $i ..."
    grep $cpg meth_CGNumber.txt >> $OUT
    let "i++"
done

echo 'DONE!'

