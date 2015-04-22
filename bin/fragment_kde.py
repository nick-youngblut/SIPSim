#!/usr/bin/env python

#--- Option parsing ---#
"""
fragment_kde: make a 2d kernel density estimate of fragment G+C and length

Usage:
  fragments [options] <fragment_table> 
  fragments -h | --help
  fragments --version

Options:
  <fragment_table>  A (pickled) table of fragment GC and lengths.
                    '-' if input from STDIN. 
  --bw=<bw>         The bandwidth scalar or function passed to
                    scipy.stats.gaussian_kde().
  -h --help         Show this screen.
  --version         Show version.
  --debug           Debug mode

Description:
  Create a 2D kernel density estimate from fragment G+C and length values
  simulated by the 'fragments' subcommand.
  
  Pickled kde objects are written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
import cPickle as pickle 

## 3rd party

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from FragGC import Frag_multiKDE


# functions
def main(args):    
    kde2d = Frag_multiKDE(args['<fragment_table>'], bandwidth=args['--bw'])
    pickle.dump(kde2d, sys.stdout)

    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    main(args)
    

        
