#!/home/jrca253/anaconda2/bin/python
import pandas as pd
import numpy as np

print "Loading GB36_to_CGNumber_dict.txt ..."
df_GB36_to_CGNumber = pd.read_table(
    "GB36_to_CGNumber_dict.txt",
    header = 0,
    sep = '\t'
)


print "Converting to dictionary ..."
GB36_to_CGNumber_dict = df_GB36_to_CGNumber.set_index('CHR:POS').T.to_dict('list')
del df_GB36_to_CGNumber


print "Loading meth_hg18.txt ..."
df = pd.read_table(
    'meth_hg18.txt',
    header = 0,
    sep = '\t'
)


print "Converting from CHR:POS to CG Number ..."
df['position'] = df['position'].map(GB36_to_CGNumber_dict)
df = df[ df['position'].notnull() ]


df.to_csv(
    'meth_GCNumber.txt',
    header = True,
    index = False,
    sep = '\t'
)

