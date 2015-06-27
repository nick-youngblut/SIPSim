#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_sampleData -- make a 'sample_data' table (phyloseq) from the OTU table

Usage:
  OTU_sampleData [options] <OTU_table_file>
  OTU_sampleData -h | --help
  OTU_sampleData --version

Options:
  <OTU_table_file>    OTU table file name.
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode

Description:
  Create a sample_data table for import into phyloseq.
  The sample_data table will consist of BD value stats for each
  sample (fraction).
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from OTU_Table import OTU_table

import numpy as np

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    
    # loading otu table
    otu_tbl = OTU_table.from_csv(Uargs['<OTU_table_file>'], sep='\t')

    # otu table in long format
    if not otu_tbl.is_long:
        otu_tbl.wide2long()

    # editing columns
    df = otu_tbl.df[['library','fraction']].drop_duplicates()
    L = lambda x: x['library'] + '__' + x['fraction']
    df['sample'] = df.apply(L, 1)
    df[['BD_min','BD_max']] = df['fraction'].str.extract('(.+)-(.+)')
    L = lambda x: round((float(x['BD_min']) + float(x['BD_max'])) / 2, 4)
    df['BD_mean'] = df.apply(L, 1)
    cols = df.columns.tolist()
    cols = ['sample'] + [x for x in cols if x != 'sample']
    df = df[cols]
    df.sort(['sample'], inplace=True)    

    # writing otu table
    df.to_csv(sys.stdout, sep='\t', index=False)                 
