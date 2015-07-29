## batteries
import sys, os
import time
import cPickle as pickle 
from functools import partial
## 3rd party
import numpy as np
import scipy.stats as stats
import dill
from pathos.multiprocessing import ProcessingPool
## Application
import Utils


# functions
def kde_intersect(x, **kwargs):
    """Wrapper for unpacking dict objects containing
    KDEs for the provided taxon.
    Args:
    taxon -- str; taxon name
    d1 -- {taxon_name:kde}
    d2 -- {taxon_name:kde}
    Returns:
    list -- [taxon_name, BD_shift]
    """
    assert len(x) >= 3, 'x must contain (taxon,kde1,kde2)'
    taxon = x[0]    

    sys.stderr.write('  Processing: {}\n'.format(taxon))
    y = _kde_intersect(x[1], x[2], **kwargs)
    assert np.isnan(y) or 0 <= y <= 1, \
        'KDE intersection is not in 0-1 range; value={}'.format(y)
    return [taxon, 1 - y]


def _kde_intersect(kde1, kde2, start=1.66, end=1.85, step=0.001):
    """Calculating the intersection of 2 KDE objects.
    np.trapz is used for integration.
    start, end, & step define the points evaluated for each KDE.
    Args:
    kde1 -- scipy kde object
    kde2 -- scipy kde object
    start -- float; start of series
    end -- float; end of series
    step -- float; step size of series
    Returns:
    float -- intersection value
    """
    if kde1 is None or kde2 is None:
        return np.NAN

    # evalution grid
    x = np.arange(start,end,step)
    # calculate intersection densities
    pmin = np.min(np.c_[kde1(x),kde2(x)], axis=1)
    # integrate areas under curves
    total = np.trapz(y=kde1(x), x=x) + \
            np.trapz(y=kde2(x), x=x)
    intersection = np.trapz(y=pmin,x=x)

    # overlap coefficient
    return 2 * intersection / float(total)


def is_kde_lib(d):
    """Is the kde object a dict of dicts {lib:{taxon:scipy_kde}}?
    Args:
    d -- {taxon_name:kde} or {libID:{taxon_name:kde}}
    Returns:
    boolean
    """
    try:
        k1 = d.keys()[0]
        d[k1].keys()
    except AttributeError:
        return False
    else:
        return True


def kde_add_lib(d):    
    """Adding top-level library ID ('NA') to dict if not present.
    Args:
    d -- {taxon_name:kde} or {libID:{taxon_name:kde}}
    Returns:
    dict -- {libID:{taxon_name:kde}}
    """
    is_lib = is_kde_lib(d)
    if is_lib:
        return d
    else:
        return {'NA':d}

        
def taxon_overlap(d1, d2):
    """Determine which taxa overlap between KDE objects.
    Warnings for taxa that do not overlap.
    Args:
    d1 -- {taxon_name:kde}
    d2 -- {taxon_name:kde}
    Returns:
    iterable -- [overlapping taxon names]
    """
    taxa1 = set(d1.keys())
    taxa2 = set(d2.keys())
    
    msg = 'WARNING: taxon "{}" not in both KDEs. It will be skipped.\n'
    for x in taxa1 ^ taxa2:
        sys.stderr.write(msg.format(x))
    
    return taxa1 & taxa2


def main(args):
    """Main function for calculating BD shift.
    """
    sys.stderr.write('Loading KDE objects...\n')
    kde1 = Utils.load_kde(args['<kde1>'])
    kde2 = Utils.load_kde(args['<kde2>'])

    # adding top-level library ID if not present
    kde1 = kde_add_lib(kde1)
    kde2 = kde_add_lib(kde2)

    sys.stderr.write('Calculating BD shifts...\n')
    print '\t'.join(['lib1','lib2','taxon','BD_shift'])
    for libID1,d1 in kde1.items():
        for libID2,d2 in kde2.items():
            msg = '  Comparing libraries; "{}", "{}"\n'
            sys.stderr.write(msg.format(libID1, libID2))

            # overlap of taxa btw libraries
            taxa = taxon_overlap(d1, d2)            

            # calculating BD shift (in parallel)
            pfunc = partial(kde_intersect, 
                            start=float(args['--start']),
                            end=float(args['--end']),
                            step=float(args['--step']))

            pool = ProcessingPool(nodes=int(args['--np']))
            if args['--debug']:
                res = map(pfunc, [(taxon, d1[taxon], d2[taxon])
                                  for taxon in taxa])
            else:
                res = pool.amap(pfunc, [(taxon, d1[taxon], d2[taxon])
                                        for taxon in taxa])
                while not res.ready():
                    time.sleep(2)
                res = res.get()        
                            
            # writing out table
            for line in res:
                print '\t'.join([libID1, libID2] + \
                                [str(x) for x in line])
                            
