#!/home/jrca253/anaconda2/bin/python
import pandas as pd

meth_file_K = "meth_K.txt"
meth_file_N = "meth_N.txt"
meth_file_T = "meth_T.txt"

out_file_K = "meth_K_cpgs_in_KNT.txt"
out_file_N = "meth_N_cpgs_in_KNT.txt"
out_file_T = "meth_T_cpgs_in_KNT.txt"


##### READ IN METH DATA #####
print "Loading " + meth_file_K + " ..."
df_K = pd.read_table(
    meth_file_K,
    sep = '\t',
    header = 0
)

print "Loading " + meth_file_N + " ..."
df_N = pd.read_table(
    meth_file_N,
    sep = '\t',
    header = 0
)

print "Loading " + meth_file_T + " ...\n\n"
df_T = pd.read_table(
    meth_file_T,
    sep = '\t',
    header = 0
)


##### KEEP ONLY CpGs PRESENT IN ALL DATA SETS #####
print "Reducing " + meth_file_K + " ..."
df_Kfinal = df_K.loc[ df_K['position'].isin(df_N['position']) & df_K['position'].isin(df_T['position']) ]
print "Reducing " + meth_file_N + " ..."
df_Nfinal = df_N.loc[ df_N['position'].isin(df_K['position']) & df_N['position'].isin(df_T['position']) ]
print "nReducing " + meth_file_T + " ..."
df_Tfinal = df_T.loc[ df_T['position'].isin(df_K['position']) & df_T['position'].isin(df_N['position']) ]

print"\n\nK:"
print df_Kfinal.head()
print"\nN:"
print df_Nfinal.head()
print"\nT:"
print df_Tfinal.head()

print "\n\nK CpGs: " + str(df_Kfinal.shape[0])
print "N CpGs: " + str(df_Nfinal.shape[0])
print "T CpGs: " + str(df_Tfinal.shape[0])


##### SAVE #### 
print "\n\nSaving " + out_file_K + " ..."
df_Kfinal.to_csv(
    out_file_K,
    sep = '\t',
    header = True,
    index = False
)

print "Saving "+ out_file_N + " ..."
df_Nfinal.to_csv(
    out_file_N,
    sep= '\t',
    header = True,
    index = False
)

print "Saving "+ out_file_T + " ..."
df_Tfinal.to_csv(
    out_file_T,
    sep= '\t',
    header = True,
    index = False
)
