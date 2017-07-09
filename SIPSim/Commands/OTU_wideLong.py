#!/usr/bin/env python

"""
OTU_wideLong: convert OTU table from wide to long or vice versa

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

  Output
  ------
  A tab-delimieted OTU table written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
from SIPSim.OTU_Table import OTU_table


def main(args=None):    
    otu_tbl = OTU_table.from_csv(args['<OTU_table_file>'], sep='\t')

    index = False
    if args['--wide']:
        otu_tbl.long2wide(values='count',index='taxon',
                          columns=['library','fraction'])
        index = True
        
    otu_tbl.to_csv(sys.stdout, sep='\t', index=index)                 


def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
