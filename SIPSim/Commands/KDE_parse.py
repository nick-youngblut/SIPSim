#!/usr/bin/env python

"""
KDE_parse: parse out KDEs for certain taxa

Usage:
  KDE_parse [options] <kde> <taxa>
  KDE_parse -h | --help
  KDE_parse --version

Options:
  <kde>         Pickled KDE object.
                ('-' if input from STDIN)
  <taxa>        List of taxa used for parsing (one name per line).
                Anything following a <tab> will be ignored.
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
import dill
## application libraries
#scriptDir = os.path.dirname(__file__)
#libDir = os.path.join(scriptDir, '../lib/')
#sys.path.append(libDir)
# application
from SIPSim import Utils
    


def load_taxa(inFile):
    """Loading a file listing taxa names.
    """
    taxa = []
    with open(inFile, 'rb') as inFH:
        for line in inFH:
            line = line.rstrip().split('\t')
            taxa.append(line[0])
    return taxa


    return kde_type


# main
def main(args=None):
    # loading taxa names
    taxa = load_taxa(args['<taxa>'])

    # loading KDEs
    KDEs = Utils.load_kde(args['<kde>'])
    
    # parsing KDEs 
    kde_type = Utils.KDE_type(KDEs)

    # parsing KDE
    if kde_type == 1:
        KDEs_p = [[t,k] for t,k in KDEs if t in taxa]
    elif kde_type == 2:
        KDEs_p = {t:k for t,k in KDEs.items() if t in taxa}
    elif kde_type == 3:
        KDEs_p = {}
        for libID,v in KDEs_p.items():
            KDEs_pp = {t:k for t,k in v.items() if t in taxa}
            KDEs_p[libID] = KDEs_pp
            KDEs_pp = None
    elif kde_type == 4:
        KDEs_p = {}
        for libID,filename in KDEs.items(): 
            KDE_bylib = Utils.load_kde(filename)
            KDE_bylib = {t:k for t,k in KDE_bylib.items() if t in taxa}
            KDEs_p[libID] = KDE_bylib
            KDE_bylib = None
    else:
        raise TypeError, 'KDE object type not recognized'

    # writing
    dill.dump(KDEs_p, sys.stdout)


def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)

