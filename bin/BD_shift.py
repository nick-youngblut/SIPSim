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
  -h --help      Show this screen.
  --version      Show version.
  --debug        Debug mode (no parallel processes)

Description:

  **Output**
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
import Utils


# functions
def kde_intersect(kde1, kde2, start=1.66, end=1.85, step=0.001):
    if kde1 is None or kde2 is None:
        return np.NAN

    # evalution grid
    x = np.arange(start,end,step)
    # calculate intersection densities
    pmin = np.min(np.c_[kde1(x),kde2(x)], axis=1)
    # integrate areas under curves
#    total = kde1.integrate_box_1d(start,end) + \
#            kde2.integrate_box_1d(start,end)
    total = np.trapz(y=kde1(x), x=x) + \
            np.trapz(y=kde2(x), x=x)
    intersection = np.trapz(y=pmin,x=x)

    # overlap coefficient
    return 2 * intersection / float(total)


def is_kde_lib(kde):
    """Is the kde object a dict of dicts {lib:{taxon:scipy_kde}}?
    """
    try:
        k1 = kde.keys()[0]
        kde[k1].keys()
    except AttributeError:
        return False
    else:
        return True

def kde_add_lib(kde):    
    is_lib = is_kde_lib(kde)
    if is_lib:
        return kde
    else:
        return {'NA':kde}

        
def taxon_overlap(d1, d2):
    taxa1 = set(d1.keys())
    taxa2 = set(d2.keys())
    
    msg = 'WARNING: taxon "{}" not in both KDEs. It will be skipped.\n'
    for x in taxa1 ^ taxa2:
        sys.stderr.write(msg.format(x))
    
    return taxa1 & taxa2


def main(args):
    sys.stderr.write('Loading kde objects...\n')
    kde1 = Utils.load_kde(args['<kde1>'])
    kde2 = Utils.load_kde(args['<kde2>'])

    kde1 = kde_add_lib(kde1)
    kde2 = kde_add_lib(kde2)

    sys.stderr.write('Calculating BD shifts...\n')
    for libID1,v1 in kde1.items():
        for libID2,v2 in kde2.items():
            taxa = taxon_overlap(v1, v2)
            for taxon_name in taxa:
#                taxon_name = 'Kosmotoga_olearia_TBF_19_5_1';                

                isect = kde_intersect(v1[taxon_name],
                                      v2[taxon_name])
                out = [libID1, libID2, taxon_name,
                       1 - isect]
                out = [str(x) for x in out]
                print '\t'.join(out)
#                sys.exit()
                

# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    main(args)
    
    

        
