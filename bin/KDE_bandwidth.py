#!/usr/bin/env python

#--- Option parsing ---#
"""
KDE_bandwidth: get the bandwidth of each KDE

Usage:
  KDE_bandwidth [options] <kde>
  KDE_bandwidth -h | --help
  KDE_bandwidth --version

Options:
  <kde>         Pickled KDE object
                ('-' if input from STDIN) 
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Get the KDE bandwidth (KDE factor) for each KDE in 
the provided kde object. 

  **Output**
  Tab-delim file: <taxon><tab><bandwidth>
"""

# import
## batteries
from docopt import docopt
import sys,os
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
    
    for taxon,kde in KDEs.items():
        kde_factor = 'NA' if kde is None else str(kde.factor)
        print '\t'.join([taxon, kde_factor])
    

        
