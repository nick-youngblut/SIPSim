#!/usr/bin/env python

#--- Option parsing ---#
"""
KDE_info: get info on KDE object files

Usage:
  KDE_info [options] <kde>
  KDE_info -h | --help
  KDE_info --version

Options:
  <kde>         Pickled KDE object
                ('-' if input from STDIN) 
  -n            Write the number of KDEs.
  -l            Write list of taxa names.
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Get some basic information on KDE object files.

  Output
  ------
  Values written to STDOUT
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
    
    # number of KDEs
    if args['-n']:
        try:
            print len(KDEs.keys())
        except AttributeError:
            print len(KDEs)
            
    # listing KDEs
    if args['-l']:
        try:            
            for x in KDEs.keys():
                print x
        except AttributeError:
            for x in KDEs:
                print x[0]

        
