#!/home/jrca253/anaconda2/bin/python
import pandas as pd

df = pd.read_table(
    'HumanMethylation450_15017482_v1-2.csv',
    header = 7,
    sep = ',',
    usecols = ['IlmnID', 'Chromosome_36', 'Coordinate_36'],
    dtype = object
)


df['CHR:POS'] = 'chr' + df['Chromosome_36'].map(str) + ':' + df['Coordinate_36'].map(str)
df = df[ ['CHR:POS', 'IlmnID'] ]


df.to_csv(
    'GB36_to_CGNumber_dict.txt',
    header = True,
    index = False,
    sep = '\t'
)

