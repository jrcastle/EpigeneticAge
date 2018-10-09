#!/home/jrca253/anaconda2/bin/python
import pandas as pd
import random
random.seed(123) 

train_fraction = 0.8
vali_fraction  = 1.0 - train_fraction

cov_filename  = "cov_K.txt"
meth_filename = "meth_K_imputed.txt"


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
meth_columns = df_meth.columns.values.tolist()

if meth_columns[1:] != cov_columns[1:]:
    print "Samples in %s and %s do not match!" % (cov_filename, meth_filename)
    print "Please resolve this discrepancy before running this code"
    exit()
#ENDIF


##### SPLIT COLUMNS RANDOMLY #####
FIRST = True
cov_train_samples  = []
cov_vali_samples   = [] 
meth_train_samples = []
meth_vali_samples  = []

it = 0
print "Splitting ..." 
for sample in cov_columns:
    
    if it % 20 == 0:
        print "Processing sample %i of %i ..." % (it, len(cov_columns))
    #ENDIF

    if FIRST:
        cov_train_samples.append( cov_columns[0] )
        cov_vali_samples.append( cov_columns[0] )
        meth_train_samples.append( meth_columns[0] )
        meth_vali_samples.append( meth_columns[0] )
        FIRST = False
        it = it + 1
        continue
    #ENDIF

    if random.random() <= train_fraction:
        cov_train_samples.append( sample )
        meth_train_samples.append( sample )
    else:
        cov_vali_samples.append( sample )
        meth_vali_samples.append( sample )
    #ENDIFELSE

    it = it + 1
#ENDFOR

df_cov_train  = df_cov[ cov_train_samples ]
df_cov_vali   = df_cov[ cov_vali_samples ]
df_meth_train = df_meth[ meth_train_samples ]
df_meth_vali  = df_meth[ meth_vali_samples ]


##### WRITE FILES #####
print "Writing " + cov_filename.split(".")[0] + "_train.txt" + " ..." 
df_cov_train.to_csv(cov_filename.split(".")[0] + "_train.txt", sep = '\t', index = False, header = True, na_rep = 'NA')

print "Writing " + cov_filename.split(".")[0] + "_vali.txt" + " ..."
df_cov_vali.to_csv(cov_filename.split(".")[0] + "_vali.txt", sep = '\t', index = False, header = True, na_rep = 'NA')

print "Writing " + meth_filename.split(".")[0] + "_train.txt" + " ..."
df_meth_train.to_csv(meth_filename.split(".")[0] + "_train.txt", sep = '\t', index = False, header = True, na_rep = 'NA')

print "Writing " + meth_filename.split(".")[0] + "_vali.txt" +" ..."
df_meth_vali.to_csv(meth_filename.split(".")[0] + "_vali.txt", sep = '\t', index = False, header = True, na_rep = 'NA')
