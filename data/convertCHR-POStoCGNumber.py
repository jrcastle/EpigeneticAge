#!/home/jrca253/anaconda2/bin/python
import pandas as pd
import numpy as np

meth_file  = "meth_T_AllCpGs.txt"
dict_file  = "CHR-POS_to_CGNumber_850k_dict.txt"
overlap21k = False

print "Loading POS_to_CGNumber_dict.txt ..."
df_CHR_POS_to_CGNumber = pd.read_table(
    dict_file,
    header = 0,
    sep = '\t'
)


print "Converting to dictionary ..."
CHR_POS_to_CGNumber_dict = df_CHR_POS_to_CGNumber.set_index('CHR:POS').T.to_dict('list')
del df_CHR_POS_to_CGNumber


print "Loading " + meth_file + " ..."
df = pd.read_table(
    meth_file,
    header = 0,
    sep = '\t'
)


print "Converting from CHR:POS to CG Number ..."
df['position'] = df['position'].map(CHR_POS_to_CGNumber_dict)
df = df[ df['position'].notnull() ]
df['position'] = pd.DataFrame( df['position'].values.tolist(), index = df['position'].index)
df = df[ df['position'].notnull() ]


if overlap21k:
    print "Checking if all 21k CpGs used by Horvath are in this data set"
    print "Missing CpGs will be added at the end with an NA for all methylation values"
    df_H = pd.read_table(
        'cpgs_21kdatMethUsed.txt',
        header = None,
        names = ['position']
    )

    df_H = df_H[ ~df_H['position'].isin(df['position']) ].copy()
    df   = pd.concat([df, df_H], sort = False)
#ENDIF

print "Saving tmp.csv ..."
df.to_csv(
    'tmp.csv',
    header = True,
    index = False,
    sep = ',',
    na_rep = "NA"
)
