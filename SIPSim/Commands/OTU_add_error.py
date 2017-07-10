#!/usr/bin/env python

"""
OTU_add_error: adding error to abundance values

Usage:
  OTU_add_error [options] <OTU_table>
  OTU_add_error -h | --help
  OTU_add_error --version

Options:
  <OTU_table>           OTU table file.
  --error_dist=<dc>     Distribution to use. 
                        (see numpy.random for possible distributions)
                        [Default: negative_binomial]
  --error_dist_p=<dp>   Distribution parameters (see numpy.random).
                        [Default: p:0.5]
  --version             Show version.
  --debug               Debug mode.
  -h --help             Show this screen.

Description:
  Adding sampling error to OTU count and rel_abund values in the OTU table.
  OTU counts are used to set the distribution mean parameter (eg., lambda or mu)
  so that each OTU count (with error) is drawn from a distinct error
  distribution determined by the OTU's original count value.

  Output
  ------
  A tab-delimited OTU table written to STDOUT
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
from SIPSim.Utils import parseKeyValueString as distParamParse
from SIPSim.OTU_Table import OTU_table
from SIPSim.Error_Dist import error_dist


def main(args=None):
    # parsing dist params
    args['--error_dist_p'] = distParamParse(args['--error_dist_p'])

    # loading OTU table 
    otu_tbl = OTU_table.from_csv(args['<OTU_table>'], sep='\t')

    # error distribution function
    e_dist = error_dist(args['--error_dist'],
                        args['--error_dist_p'])

    # drawing error from OTU counts
    f = lambda x : e_dist.sample(1, x)[0]
    otu_tbl.apply_each_taxon(f, 'count')

    # Setting relative abundance
    otu_tbl.add_rel_abund('count', 'rel_abund')
    
    # writing out table
    otu_tbl.to_csv(sys.stdout, sep='\t', index=False)


# main
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
    

        
