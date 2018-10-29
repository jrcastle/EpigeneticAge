#!/home/jrca253/anaconda2/bin/python
import pandas as pd

df = pd.read_table(
    'cov_K.txt',
    header = 0,
    sep = '\t'
)

df.set_index('cov',inplace=True)

df = df.T
df['Female'] = 1
df['Tissue'] = "Breast"
df = df[ ['Age', 'Female', 'Tissue'] ]
df.index.name = "Sample"

df.to_csv(
    'sample_annotation.csv',
    header = True,
    index = True,
    sep = ','
)

