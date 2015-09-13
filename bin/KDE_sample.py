#!/usr/bin/env python

#--- Option parsing ---#
"""
KDE_sample: sample from each KDE and write a table of values

Usage:
  KDE_sample [options] <kde>
  KDE_sample -h | --help
  KDE_sample --version

Options:
  <kde>         Pickled KDE object
                ('-' if input from STDIN) 
  -n=<n>        Sample size.
                [Default: 10000]
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Sample values from each KDE in the pickled KDE object
  and produce a table of values.

  Output
  ------
  Tab-delim file: <taxon><tab><value>
"""

# import
## batteries
from docopt import docopt
import sys,os
## 3rd party
import pandas as pd
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
# application
import Utils
    
    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    KDEs = Utils.load_kde(args['<kde>'])
    
    n = int(args['-n'])

    vals = {taxon:kde.resample(n)[0,] for taxon,kde in KDEs.items() \
            if kde is not None}
    tbl = pd.DataFrame(vals)
    tbl.to_csv(sys.stdout, sep='\t')
        

        
