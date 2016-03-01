#!/usr/bin/env python

#--- Option parsing ---#
"""
KDE_selectTaxa: Select a set of taxa to be incorporat

Usage:
  KDE_selectTaxa [options] <kde>
  KDE_selectTaxa -h | --help
  KDE_selectTaxa --version

Options:
  <kde>         Pickled KDE object
                ('-' if input from STDIN) 
  -s=<s>        Subsample from taxa names.
  -r            Subsample with replacement.
  -a            Select taxa if if they are missing 
                in some libraries (gradients).
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Select taxa names from the KDE object. A subsample of taxa can be selected.
  This is useful for selecting consistent set of taxa to be incorporators 
  for replicate simulations (ie., `taxa` option in isotope_incorp subcommand).
  
  Note: if -s > total_number_taxa, then subsampling will be done with 
  replacement.

  Output
  ------
  Values written to STDOUT
"""

# import
## batteries
from docopt import docopt
import sys,os
import pandas as pd
import numpy as np
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
# application
import Utils
    

# functions
def get_taxa_OLD(KDEs, all_taxa=None):
    if all_taxa is not None:
        try: 
            taxa = [k for k,v in KDEs.items()]
        except AttributeError:
            taxa = [k[0] for k in KDEs if k[1]]
    else:
        try: 
            taxa = [k for k,v in KDEs.items() if v is not None]
        except AttributeError:
            taxa = [k[0] for k in KDEs if k[1] is not None]
    return taxa


def get_taxa(KDEs, kde_type, all_taxa=None):
    """Getting taxa names from KDE object.
    """
    # parsing KDE
    if kde_type == 1: 
        if all_taxa is not None:
            taxa = [k[0] for k in KDEs]            
        else:
            taxa = [k[0] for k in KDEs if k[1] is not None]
    elif kde_type == 2:
        if all_taxa is not None:
            taxa = KDEs.keys()        
        else:
            taxa = [k for k,v in KDEs.items() if v is not None]
    elif kde_type == 3:
        taxa = []
        for libID,v in KDEs.items():
            if all_taxa is not None:
                taxa += v.keys()  
            else:
                taxa += [k for k,vv in v.items() if v is not None]
        taxa = list(set(taxa))            
    else:
        raise TypeError, 'KDE object type not recognized'

    return taxa


# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    if args['-r'] is not None:
        args['-r'] = True
    else:
        args['-r'] = False

    # load KDEs object
    sys.stderr.write('Loading KDEs...\n')
    KDEs = Utils.load_kde(args['<kde>'])

    # KDE type
    kde_type = Utils.KDE_type(KDEs)            
        
    # all taxa names in KDEs
    taxa = get_taxa(KDEs, kde_type, args['-a'])
    ntaxa = len(taxa)

    # subsampling (if needed)
    if args['-s'] is not None:
        nsub = int(args['-s'])
        if nsub > ntaxa:
            args['-r'] = True
        taxa = np.random.choice(taxa, size=nsub, replace=args['-r'])
        msg = 'Subsampled {} taxa'
        sys.stderr.write(msg.format(nsub) + '\n')
    
    # writing
    for x in taxa:
        print x
