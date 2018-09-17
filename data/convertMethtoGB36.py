#!/home/jrca253/anaconda2/bin/python
import pandas as pd
import numpy as np
from pyliftover import LiftOver
lo = LiftOver('hg19', 'hg18') # GRCh37/hg19 -> NCBI36/hg18

def convert( locus ):
    
    CHR = str(locus).split(':')[0]
    POS = int( str(locus).split(':')[1] )
    POS = POS - 1
    CONV_LIST = lo.convert_coordinate(CHR, POS)
    CONV = np.nan
    if len(CONV_LIST) == 1:
        CONV = str( CONV_LIST[0][0] ) + ':' + str( int(CONV_LIST[0][1]) + 1 )
    #ENDIF
    return CONV
#END convert


df = pd.read_table(
    'meth.txt',
    header = 0,
    sep = '\t'
)

df['position'] = df['position'].map(convert)
df = df[ df['position'].notnull() ]


df.to_csv(
    'meth_hg18.txt',
    header = True,
    index = False,
    sep = '\t'
)

