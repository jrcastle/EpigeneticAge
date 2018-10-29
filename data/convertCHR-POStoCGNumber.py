#!/home/jrca253/anaconda2/bin/python
import pandas as pd
import numpy as np

meth_file = "meth_K.txt"

print "Loading POS_to_CGNumber_dict.txt ..."
df_CHR_POS_to_CGNumber = pd.read_table(
    "CHR-POS_to_CGNumber_dict.txt",
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


print "Checking to if all Horvath Clock CpGs are in this data set"
print "Missing clock CpGs will be added at the end with an NA for all methylation values"
df_H = pd.read_table(
    'HorvathClockCpGs.txt',
    header = None,
    names = ['position']
)

#df_H = df_H[ ~df_H['position'].isin(df['position']) ].copy()
#df   = pd.concat([df, df_H], sort = False)


print "Saving meth_CGNumber.txt"
df.to_csv(
    'meth_CGNumber.txt',
    header = True,
    index = False,
    sep = ',',
    na_rep = "NA"
)

