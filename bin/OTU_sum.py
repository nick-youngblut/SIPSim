#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_sum -- Sum OTU counts (by group)

Usage:
  OTU_summary [options] <OTU_table_file>
  OTU_summary -h | --help
  OTU_summary --version

Options:
  <OTU_table_file>  Name of file produced by OTU_sim subcommand.
  --groupby=<gb>    Comma-delimited list of columns to group by.
                    [default: library,fractions]
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
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from OTU_Table import OTU_table


# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')

    groupby = Uargs['--groupby'].replace(' ','').split(',')
    
    otu_tbl = OTU_table.from_csv(Uargs['<OTU_table_file>'], sep='\t')

    try:
        otu_tbl.df.groupby(groupby).sum().to_csv(sys.stdout, sep='\t')
    except KeyError:
        sys.exit('ERROR: --groupby values (at least some) are not in the OTU table')
