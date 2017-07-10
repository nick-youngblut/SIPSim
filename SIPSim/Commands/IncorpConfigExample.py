#!/usr/bin/env python

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
  --n_reps=<nr>            Number of replicates per control/treatment community
                           [Default: 1]
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

  OUTPUT
  ------
  Config file written to STDOUT
"""

# import
## batteries
from docopt import docopt
import sys,os
import dill

## application libraries
from SIPSim import Config
from SIPSim.Config import ExampleConfig


    
def main(args=None):
    """Main entry point for command
    """
    n_reps = int(args['--n_reps'])    
    for i in xrange(n_reps):
        # creating basic config object
        basicConfig = Config.get_basicConfig()
        cfg = ExampleConfig(basicConfig)

        # editing config object
        if args['--percTaxa']:
            cfg.set_percTaxa(args['--percTaxa'])
        if args['--percIncorpUnif']:
            cfg.set_percIncorpUnif(args['--percIncorpUnif'])
        if args['--percIncorpMean']:
            if args['--percIncorpSD'] is None:
                args['--percIncorpSD'] = 0.1
            cfg.set_percIncorpNorm(args['--percIncorpMean'],
                                   args['--percIncorpSD'])

        # renaming keys & adding to cfgs collection
        for k in cfg.keys():
            newname = str(int(k) + i * 2)
            cfg.rename(k, newname)

        # writing config file
        print '\n'.join(cfg.write())

    
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)

