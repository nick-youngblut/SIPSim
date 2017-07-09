#!/usr/bin/env python

"""
fragment_KDE: make a 2d kernel density estimate of fragment 
              buoyant density and length

Usage:
  fragment_KDE [options] <fragment_table> 
  fragment_KDE -h | --help
  fragment_KDE --version

Options:
  <fragment_table>  A (pickled) table of fragment GC and lengths.
                    '-' if input from STDIN. 
  --bw=<bw>         The bandwidth scalar or function passed to
                    scipy.stats.gaussian_kde(); 'bw_method'
  -h --help         Show this screen.
  --version         Show version.
  --debug           Debug mode

Description:
  Create a 2D kernel density estimate from fragment G+C and length values
  simulated by the 'fragments' subcommand.

  Output
  ------  
  Pickled kde object written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
## 3rd party
import dill
# application
from SIPSim import Fragments as Frags


# functions
def main(args):    
    frag_tbl = Frags.load_frags(args['<fragment_table>'])    
    frag_kde = Frags.fit_kde(frag_tbl, bw_method=args['--bw'])
    dill.dump(frag_kde, sys.stdout)

    
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)


        
