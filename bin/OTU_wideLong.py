#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_wideLong -- convert OTU table from wide to long or vice versa

Usage:
  OTU_wideLong [options] <OTU_table_file>
  OTU_wideLong -h | --help
  OTU_wideLong --version

Options:
  <OTU_table_file>    OTU table file.
  -w --wide           Convert table to 'wide' format.
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode

Description:
  Re-format an OTU table (file produced from OTU_sim subcommand)
  to be either 'wide' (standard OTU table for QIIME, mothur, etc) or
  'long' format (easier for plotting and other functions).

  By default, the output table is in 'long' format.

  For the 'wide' format, column names are formatted as '[library]__[fraction]'
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


# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    
    otu_tbl = OTU_table.from_csv(Uargs['<OTU_table_file>'], sep='\t')

    index = False
    if Uargs['--wide']:
        otu_tbl.long2wide(values='count',index='taxon',
                          columns=['library','fraction'])
        index = True
        
    otu_tbl.to_csv(sys.stdout, sep='\t', index=index)                 
