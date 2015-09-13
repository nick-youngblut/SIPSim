#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_PCR: simulate PCR of gradient fraction DNA samples

Usage:
  OTU_PCR [options] <OTU_table>
  OTU_PCR -h | --help
  OTU_PCR --version

Options:
  <OTU_table>           OTU table file.
  --error_dist=<dc>     Distribution to use. 
                        (see Description for possible distributions)
                        [Default: neg_binom]
  --error_dist_p=<dp>   Distribution parameters.
                        (see Description for possible distributions)
                        [Default: alpha:0.5]
  --version             Show version.
  --debug               Debug mode.
  -h --help             Show this screen.

Description:
  Adding sampling error to OTU count and rel_abund values in the OTU table.
  OTU counts are used to set the distribution mean parameter (eg., lambda or mu)
  so that each OTU count (with error) is drawn from a distinct error
  distribution determined by the OTU's original count value.

  Error distributions
  -------------------
    neg_binom : negative binomial distribution (params: alpha)
        alpha param sets distribution shape.
    others : all other distributions and parameters are obtained from numpy.random

  Output
  ------
  A tab-delimited OTU table written to STDOUT
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from Utils import parseKeyValueString as distParamParse
from OTU_Table import OTU_table
from Error_Dist import error_dist


def main(Uargs):
    # parsing dist params
    Uargs['--error_dist_p'] = distParamParse(Uargs['--error_dist_p'])

    # loading OTU table 
    otu_tbl = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')

    # error distribution function
    e_dist = error_dist(Uargs['--error_dist'],
                        Uargs['--error_dist_p'])

    # drawing error from OTU counts
    f = lambda x : e_dist.sample(1, x)[0]
    otu_tbl.apply_each_taxon(f, 'count')

    # Setting relative abundance
    otu_tbl.add_rel_abund('count', 'rel_abund')
    
    # writing out table
    otu_tbl.to_csv(sys.stdout, sep='\t', index=False)


# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)

    

        
