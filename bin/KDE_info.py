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
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
# application
import Utils
    

# functions
def KDE_dataset_stats(kde, ID):
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
        print '\t'.join([str(x) for x in stats[i,]])
    


# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            

    # load KDEs object
    sys.stderr.write('Loading KDEs...\n')
    KDEs = Utils.load_kde(args['<kde>'])
    
    # number of KDEs
    if args['-n']:
        try:
            print len(KDEs.keys())
        except AttributeError:
            print len(KDEs)

    # header
    if args['-s']:
        print '\t'.join(['taxon_ID', 'KDE_ID', 'min', 'percentile_5', 
                         'percentile_25', 'mean', 'median', 'percentile_75', 
                         'percentile_95', 'max', 'stdev'])
            
    # listing KDEs 
    if args['-t'] or args['-s']:
        try:            
            for x in KDEs.keys():
                if args['-s']:
                    KDE_dataset_stats(KDEs[x], x)
                else:
                    print x
                
        except AttributeError:
            for x in KDEs:
                if args['-s']:
                    KDE_dataset_stats(x[1], x[0])
                else:
                    print x[0]
                            
