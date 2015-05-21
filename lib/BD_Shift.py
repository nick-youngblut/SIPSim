## batteries
import sys, os
import cPickle as pickle 
from functools import partial
## 3rd party
import parmap
import numpy as np
import scipy.stats as stats
## Application
import Utils


# functions
def kde_intersect(taxon, d1, d2, **kwargs):
    """Wrapper for unpacking dict objects containing
    KDEs for the provided taxon.
    Args:
    taxon -- str; taxon name
    d1 -- {taxon_name:kde}
    d2 -- {taxon_name:kde}
    Returns:
    list -- [taxon_name, BD_shift]
    """
    x = _kde_intersect(d1[taxon], d2[taxon], **kwargs)
    assert np.isnan(x) or 0 <= x <= 1, \
        'KDE intersection is not in 0-1 range; value={}'.format(x)
    return [taxon, 1 - x]


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
#    total = kde1.integrate_box_1d(start,end) + \
#            kde2.integrate_box_1d(start,end)
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
            pfunc = partial(kde_intersect, d1=d1, d2=d2,
                            start=float(args['--start']),
                            end=float(args['--end']),
                            step=float(args['--step']))
            res = parmap.map(pfunc, taxa, 
                             processes=int(args['--np']),
                             chunksize=int(args['--cs']),
                             parallel=not args['--debug'])
            
            # writing out table
            for line in res:
                print '\t'.join([libID1, libID2] + \
                                [str(x) for x in line])
                            
