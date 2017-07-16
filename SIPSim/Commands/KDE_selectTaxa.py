#!/usr/bin/env python

"""
KDE_selectTaxa: Select a set of taxa to be incorporators

Usage:
  KDE_selectTaxa [options] <kde>
  KDE_selectTaxa -h | --help
  KDE_selectTaxa --version

Options:
  <kde>         Pickled KDE object
                ('-' if input from STDIN) 
  -s=<s>        Subsample `s` number of taxa.
  -f=<f>        Subsample `f` fraction of taxa. 
  -p=<p>        Subsample `p` percentage of taxa.
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
  
  If -s > total_number_taxa or -f is > 1 or -p is > 100, then subsampling will 
  be done with replacement.

  -s takes precedent over -f & -p. 

  Output
  ------
  Taxon names written to STDOUT (1 per line)
"""

# import
## batteries
from docopt import docopt
import sys,os
import pandas as pd
import numpy as np
## application libraries
from SIPSim import Utils
    

# functions
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
    elif kde_type == 4:
        taxa = []
        for libID,filename in KDEs.items(): 
            KDE_bylib = Utils.load_kde(filename)
            if all_taxa is not None:
                taxa += KDE_bylib.keys()        
            else:
                taxa += [k for k,v in KDE_bylib.items() if v is not None]
    else:
        raise TypeError, 'KDE object type not recognized'

    return taxa


def main(args=None):
    # subsample with replacement arg
    #if args['-r'] is not None:
    #    args['-r'] = True
    #else:
    #    args['-r'] = False

    # load KDEs object
    sys.stderr.write('Loading KDEs...\n')
    KDEs = Utils.load_kde(args['<kde>'])

    # KDE type
    kde_type = Utils.KDE_type(KDEs)            
        
    # all taxa names in KDEs
    taxa = get_taxa(KDEs, kde_type, args['-a'])
    ntaxa = len(taxa)

    # subsampling (if needed)
    ## number to subsample
    nsub = None
    if args['-s'] is not None:
        nsub = int(args['-s'])
    elif args['-f'] is not None:
        nsub = float(args['-f']) * ntaxa
        nsub = int(nsub)
    elif args['-p'] is not None:
        nsub = float(args['-p']) / 100 * ntaxa
        nsub = int(nsub)        
    ## subsampling
    if nsub is not None:
        if nsub > ntaxa:
            args['-r'] = True
            msg = 'WARNING: nsub > ntaxa, sub-sampling with replacement!'
            sys.stderr.write(msg + '\n')
        taxa = np.random.choice(taxa, size=nsub, replace=args['-r'])
        msg = 'Subsampled {} taxa'
        sys.stderr.write(msg.format(nsub) + '\n')
    
    # writing to STDOUT
    for x in taxa:
        print x

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
