#!/usr/bin/env python

#--- Option parsing ---#
"""
incorpConfigExample: create example isotope incorporation config file

Usage:
  incorpConfigExample [options]
  incorpConfigExample -h | --help
  incorpConfigExample --version

Options:
  --percTaxa=<pt>          Max percent taxa with any incorporation
  --percIncorpUnif=<piu>   Percent incorporation (uniform distribution)
  --percIncorpMean=<pim>   Mean percent incorp (normal distribution)
  --percIncorpSD=<pis>     Stdev percent incorp (normal distribution)
  -h --help                Show this screen.
  --version                Show version.
  --debug                  Debug mode

Description:
  The example config will have 2 communities: a 'control'
  (no isotope) incorporation and a 'treatment' some possible
  incorporation. 
  The options will modify some basic parameters of how much
  isotope is incorporated in the treatment community.

  See the isoIncorp subcommand documentation for more info
  on isotope incorporation config files.

  OUTPUT:
    Config file written to STDOUT
"""

# import
## batteries
from docopt import docopt
import sys,os
import dill

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import Config
from Config import ExampleConfig


    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')    
    
    # creating basic config object
    basicConfig = Config.get_basicConfig()
    cfg = ExampleConfig(basicConfig)

    # editing config object
    if args['--percTaxa']:
        cfg.set_percTaxa(args['--percTaxa'])
    if args['--percIncorpUnif']:
        cfg.set_percIncorpUnif(args['--percIncorpUnif'])
    if args['--percIncorpMean'] or args['--percIncorpSD']:
        if args['--percIncorpMean'] is None:
            args['--percIncorpMean'] = 90
        if args['--percIncorpSD'] is None:
            args['--percIncorpSD'] = 1
        cfg.set_percIncorpNorm(args['--percIncorpMean'],
                               args['--percIncorpSD'])

    # writing config file
    print '\n'.join(cfg.write())

    

