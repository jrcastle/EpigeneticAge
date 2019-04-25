#!/home/jrca253/anaconda2/bin/python
import pandas as pd

df = pd.read_table(
    'MethylationEPIC_v-1-0_B4.csv',
    header = 7,
    sep = ',',
    usecols = ['IlmnID', 'CHR', 'MAPINFO'],
    dtype = object
)


df['CHR:POS'] = 'chr' + df['CHR'].map(str) + ':' + df['MAPINFO'].map(str)
df = df[ ['CHR:POS', 'IlmnID'] ]

df.drop( df.loc[ df['CHR:POS'] == "chrMULTI:MULTI" ].index, inplace = True )
df.drop( df.loc[ df['CHR:POS'] == "chrnan:nan" ].index, inplace = True )


df.to_csv(
    'CHR-POS_to_CGNumber_850k_dict.txt',
    header = True,
    index = False,
    sep = '\t'
)

