#!/usr/bin/env python

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

  Output
  ------
  Tab-delim file: <taxon><tab><bandwidth>
"""

# import
## batteries
from docopt import docopt
import sys,os
## application libraries
from SIPSim import Utils
    

def kde_factor(kde):
    kde_factor = 'NA' if kde is None else str(kde.factor)
    return kde_factor 

def main(args=None):
    KDEs = Utils.load_kde(args['<kde>'])

    # KDE object type
    kde_type = Utils.KDE_type(KDEs)

    # info for KDEs
    print '\t'.join(['libID', 'taxon', 'bandwidth'])
    if kde_type == 1: 
        for t,k in KDEs:
            print '\t'.join(['1', t, kde_factor(k)])
    elif kde_type == 2:
        for t,k in KDEs.items():
            print '\t'.join(['1', t, kde_factor(k)])
    elif kde_type == 3:
        for libID,x in KDEs.items():
            for t,k in x.items():
                print '\t'.join([libID, t, kde_factor(k)])
    elif kde_type == 4:
        for libID,filename in KDEs.items(): 
            KDE_bylib = Utils.load_kde(filename)
            for t,k in KDE_bylib.items():
                print '\t'.join([libID, t, kde_factor(k)])        
    else:
        raise TypeError, 'KDE object type not recognized'

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)

