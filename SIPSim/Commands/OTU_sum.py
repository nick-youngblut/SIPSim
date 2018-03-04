#!/usr/bin/env python

"""
OTU_sum: Sum OTU counts (by group)

Usage:
  OTU_sum [options] <OTU_table_file>
  OTU_sum -h | --help
  OTU_sum --version

Options:
  <OTU_table_file>  OTU table file.
  --groupby=<gb>    Comma-delimited list of columns to group by.
                    [default: library,fraction]
  --version         Show version.
  --debug           Debug mode

Description:
  Sum OTU counts. The OTU table can be grouped by
  other columns to produce aggregate sums.
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
from SIPSim.OTU_Table import OTU_table


def main(args=None):
    groupby = args['--groupby'].replace(' ','').split(',')
    
    otu_tbl = OTU_table.from_csv(args['<OTU_table_file>'], sep='\t')

    try:
        otu_tbl.df.groupby(groupby).sum().to_csv(sys.stdout, sep='\t')
    except KeyError:
        sys.exit('ERROR: --groupby values (at least some) are not in the OTU table')


def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
          
