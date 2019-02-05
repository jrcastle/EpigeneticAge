#!/home/jrca253/anaconda2/bin/python
import pandas as pd
import numpy as np

meth_file  = "HorvathClockCpGs_CGNumber.txt"
dict_file  = "CHR-POS_to_CGNumber_dict.txt"
#dict_file  = "CHR-POS_to_CGNumber_850k_dict.txt"
overlap21k = False

print "Loading POS_to_CGNumber_dict.txt ..."
df_CHR_POS_to_CGNumber = pd.read_table(
    dict_file,
    header = 0,
    sep = '\t'
)


print "Converting to dictionary ..."
CHR_POS_to_CGNumber_dict = df_CHR_POS_to_CGNumber.set_index('IlmnID').T.to_dict('list')
del df_CHR_POS_to_CGNumber


print "Loading " + meth_file + " ..."
df = pd.read_table(
    meth_file,
    header = None,
    names = ["position"]
)


print "Converting from CG Number to CHR:POS ..."
df['position'] = df['position'].map(CHR_POS_to_CGNumber_dict)
df = df[ df['position'].notnull() ]
df['position'] = pd.DataFrame( df['position'].values.tolist(), index = df['position'].index)
df = df[ df['position'].notnull() ]


print "Saving tmp.csv ..."
df.to_csv(
    'tmp.csv',
    header = True,
    index = False,
    sep = ',',
    na_rep = "NA"
)
