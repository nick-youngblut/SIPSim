#!/usr/bin/env python

#--- Option parsing ---#
"""
fragment_kde: make a 2d kernel density estimate of fragment 
              buoyant density and length

Usage:
  fragments [options] <fragment_table> 
  fragments -h | --help
  fragments --version

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

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import Fragments as Frags


# functions
def main(args):    
    frag_tbl = Frags.load_frags(args['<fragment_table>'])    
    frag_kde = Frags.fit_kde(frag_tbl, bw_method=args['--bw'])
    dill.dump(frag_kde, sys.stdout)

    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    main(args)
    

        
