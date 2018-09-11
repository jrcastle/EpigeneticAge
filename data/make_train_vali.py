#!/home/jrca253/anaconda2/bin/python
import pandas as pd
import random

train_fraction = 0.8
vali_fraction  = 1.0 - train_fraction

cov_filename  = "cov_noNA.txt"
meth_filename = "meth_noNA_addCov.txt"
suffix        = "_noNA_addCov.txt"


##### LOAD DATA #####
print "Loading " + cov_filename + " ..."
df_cov = pd.read_table(
    cov_filename,
    header = 0,
    sep = '\t'
)
cov_columns = df_cov.columns.values.tolist()

print "Loading " + meth_filename + " ..."
df_meth = pd.read_table(
    meth_filename,
    header = 0,
    sep = '\t'
)
meth_columns = []
meth_columns.append( df_meth.columns.values.tolist()[0] )


##### RENAME METH COLUMNS TO MATCH COV COLUMNS
for i in range(1, len(cov_columns)):
    meth_columns.append( cov_columns[i] )
#ENDFOR

df_meth.columns = meth_columns


##### SPLIT COLUMNS RANDOMLY #####
FIRST = True
cov_train_samples  = []
cov_vali_samples   = [] 
meth_train_samples = []
meth_vali_samples  = []


print "Splitting ..." 
for sample in cov_columns:

    if FIRST:
        cov_train_samples.append( cov_columns[0] )
        cov_vali_samples.append( cov_columns[0] )
        meth_train_samples.append( meth_columns[0] )
        meth_vali_samples.append( meth_columns[0] )
        FIRST = False
        continue
    #ENDIF

    if random.random() <= train_fraction:
        cov_train_samples.append( sample )
        meth_train_samples.append( sample )
    else:
        cov_vali_samples.append( sample )
        meth_vali_samples.append( sample )
    #ENDIFELSE

#ENDFOR

df_cov_train  = df_cov[ cov_train_samples ]
df_cov_vali   = df_cov[ cov_vali_samples ]
df_meth_train = df_meth[ meth_train_samples ]
df_meth_vali  = df_meth[ meth_vali_samples ]


##### WRITE FILES #####
df_cov_train.to_csv("cov_train" + suffix,   sep = '\t', index = False, header = True, na_rep = 'NA')
df_cov_vali.to_csv("cov_vali" + suffix,     sep = '\t', index = False, header = True, na_rep = 'NA')
df_meth_train.to_csv("meth_train" + suffix, sep = '\t', index = False, header = True, na_rep = 'NA')
df_meth_vali.to_csv("meth_vali" + suffix,   sep = '\t', index = False, header = True, na_rep = 'NA')
