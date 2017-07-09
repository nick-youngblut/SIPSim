#!/usr/bin/env python

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
  -t            Write list of taxa names.
  -s            Write summary stats on the dataset associated with each KDE.
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
import pandas as pd
import numpy as np
## application libraries
#scriptDir = os.path.dirname(__file__)
#libDir = os.path.join(scriptDir, '../lib/')
#sys.path.append(libDir)
# application
from SIPSim import Utils
    

# functions
def KDE_dataset_stats(kde, ID, libID='1'):
    """summary stats on kde.dataset attribute
    """
    try:
        dataset_dims = kde.dataset.shape
    except AttributeError:
        dataset_dims = None

    if dataset_dims:
        # stats
        ids = [ID] * dataset_dims[0]
        KDE_ids = np.arange(dataset_dims[0]) + 1
        mins = np.min(kde.dataset, axis=1)
        perc5s = np.percentile(kde.dataset, 5, axis=1)
        q1s = np.percentile(kde.dataset, 25, axis=1)
        means = np.mean(kde.dataset, axis=1)
        medians = np.median(kde.dataset, axis=1)
        q3s = np.percentile(kde.dataset, 75, axis=1)
        perc95s = np.percentile(kde.dataset, 95, axis=1)
        maxs = np.max(kde.dataset, axis=1)
        stds = np.std(kde.dataset, axis=1)

        # combining objects
        if dataset_dims > 1:
            stats = np.vstack((ids, KDE_ids, mins, perc5s, q1s, means, medians, 
                               q3s, perc95s, maxs, stds))
            stats = np.transpose(stats)
        else:
            stats = np.concatenate(ids, KDE_ids, mins, perc5s, q1s, means, 
                                   medians, q3s, perc95s, maxs, stds)
    else:
        # No KDE: nan
        stats = np.array([[ID] + [np.nan] * 10])

    # output
    for i in xrange(stats.shape[0]):
        print '\t'.join([libID] + [str(x) for x in stats[i,]])
    


def main(args=None):

    # load KDEs object
    sys.stderr.write('Loading KDEs...\n')
    KDEs = Utils.load_kde(args['<kde>'])
    
    # header
    if args['-s']:
        print '\t'.join(['lib_ID', 'taxon_ID', 'KDE_ID', 'min', 'percentile_5', 
                         'percentile_25', 'mean', 'median', 'percentile_75', 
                         'percentile_95', 'max', 'stdev'])


    # KDE type
    kde_type = Utils.KDE_type(KDEs)

    # parsing KDE
    if kde_type == 1: 
        if args['-n']:
            print len(KDEs)
            sys.exit()
        for x in KDEs:   
            if args['-s']:
                KDE_dataset_stats(x[1], x[0])
            else:
                print x[0]        
    elif kde_type == 2:
        if args['-n']:
            print len(KDEs.keys())
            sys.exit()
        for x,y in KDEs.items():
            if args['-s']:
                KDE_dataset_stats(y, x)
            else:
                print x
    elif kde_type == 3:
        if args['-n']:
            print '\t'.join(['library', 'N'])
            for libID,v in KDEs.items():
                print '\t'.join([str(x) for x in [libID, len(v.keys())]])
            sys.exit()
        for x,y in KDEs.items():
            for xx,yy in y.items(): 
                if args['-s']:
                    KDE_dataset_stats(yy, xx, libID=x)
                else:
                    print '\t'.join([x,xx])        
    elif kde_type == 4:
        for libID,filename in KDEs.items(): 
            KDE_bylib = Utils.load_kde(filename)
            if args['-n']:
                print len(KDE_bylib.keys())
                sys.exit()
            for x,y in KDE_bylib.items():
                if args['-s']:
                    KDE_dataset_stats(y, x)
                else:
                    print x
    else:
        raise TypeError, 'KDE object type not recognized'

    
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
            
