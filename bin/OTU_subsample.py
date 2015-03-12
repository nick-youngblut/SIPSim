#!/usr/bin/env python

"""
OTU_subsample: subsample without replacement from an OTU table

Usage:
  OTU_subsample [options] <OTU_table_file>
  OTU_subsample -h | --help
  OTU_subsample --version

Options:
  <OTU_table_file>    Name of file produced by OTU_table subcommand.
  --n_dist=<nd>       Distribution used to select number of samples per community.
                      See numpy.random for possible distributions.
                      For the same number of samples per community, use 'uniform'.
                      [default: uniform]
  --n_dist_p=<dp>     Parameters for the distribution used to select the number. 
                      of samples per community.
                      Input format: 'key1:value1,key2:value2,keyN:valueN'
                      For the same number of samples per community (say, 1000):
                        use: 'low:1000,high:1000'.
                      [default: low:100,high:100]
  --samp_min          Use the minimum number of taxa in any community as the number
                      subsampled for all communities.
                      Overrides 'n_dist' and 'n_dist_p' 
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode

Description:
  Subsample from an OTU table created by the OTU_table subcommand.

  This is used to simulate the actual sequencing of DNA/RNA from each gradient fraction.
  
  n_dist and n_dist_p options can be used to simulate uneven numbers of sequences
  (here, signified as OTUs) per sample. 
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from OTU_table import OTU_table
#from SimComms import SimComms
import Utils


def main(Uargs):
    # dist params as dict
    Uargs['--n_dist_p'] = Utils.parseKeyValueString(Uargs['--n_dist_p'])

    
    otu_tbl = OTU_table.from_csv(Uargs['<OTU_table_file>'], sep='\t')

    print otu_tbl; sys.exit()
    

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)
    

        