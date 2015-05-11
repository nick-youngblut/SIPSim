#!/usr/bin/env python

#--- Option parsing ---#
"""
BD_shift: Determine the shift in BD based on KDE overlap

Usage:
  fragments [options] <kde1> <kde2>
  fragments -h | --help
  fragments --version

Options:
  <kde1>         KDE object
  <kde2>         KDE object 
  --np=<np>      Number of parallel processes.
                 [default: 1]
  --cs=<cs>      Chunksize for each process (number of taxa).
                 [default: 1]
  --start=<s>    Start of series for integration.
                 [default: 1.66]
  --end=<e>      End of series for integration.
                 [default: 1.85]
  --step=<x>     Step size of series for integration.
                 [default: 0.001]
  --version      Show version.
  --debug        Debug mode (no parallel processes)
  -h --help      Show this screen.

Description:
  Determine the shift in BD value distribution between 2
  KDEs of BD values.
  The BD shift is calculated as 1 - KDE_intersection, 
  where KDE_intersection is calculated from 
  the integral of intersection densities at each point
  in the specified value series (start,end,step).
 
  **Output**
  A tab-delimited table of BD shift values is written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
import cPickle as pickle 
from functools import partial
## 3rd party
import parmap
import numpy as np
import scipy.stats as stats

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import BD_Shift


# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    BD_Shift.main(args)
    
    

        
