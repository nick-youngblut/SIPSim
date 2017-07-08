## batteries
import sys, os
import time
import cPickle as pickle 
from functools import partial
## 3rd party
from docopt import docopt
import numpy as np
import scipy.stats as stats
import dill
from pathos.multiprocessing import ProcessingPool
## Application
import Utils


# functions
def kde_intersect(x, **kwargs):
    """Wrapper for unpacking dict objects containing KDEs for the
    provided taxon.
 
    Parameters
    ----------
    x : list
        [taxon name, kde1, kde2]
    kwargs : optional
        passed to `_kde_intersect`

    Returns
    -------
    list
       [taxon_name, BD_shift]

    Notes
    -----
    The BD shift values range from 0 to 1 (1 = no BD distribution overlap).
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
    (start, end, & step) define the points evaluated for each KDE.
    
    Parameters
    ----------
    kde1, kde2 : scipy kde object
    start : float, optional
        start of BD value series
    end :  float, optional
        end of BD value series
    step : float, optional
        step size of series

    Returns
    -------
    intersect_value : float
        overal coefficient
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

    Parameters
    ----------
    d : dict
        {taxon_name:kde} or {libID:{taxon_name:kde}}
    Returns
    -------
    b : boolean
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

    Parameters
    ----------
    d : dict
        {taxon_name:kde} or {libID:{taxon_name:kde}}

    Returns
    -------
    d : dict 
        {libID:{taxon_name:kde}}
    """
    is_lib = is_kde_lib(d)
    if is_lib:
        return d
    else:
        return {'NA':d}

        
def taxon_overlap(d1, d2):
    """Determine which taxa overlap between KDE objects.
    Warnings for taxa that do not overlap.

    Parameters
    ----------
    d1, d2 : dict
        {taxon_name:kde}

    Returns
    -------
    taxon_names : iterable
        set of overlapping taxon names
    """
    taxa1 = set(d1.keys())
    taxa2 = set(d2.keys())
    
    msg = 'WARNING: taxon "{}" not in both KDEs. It will be skipped.\n'
    for x in taxa1 ^ taxa2:
        sys.stderr.write(msg.format(x))
    
    return taxa1 & taxa2


def main(args):
    """Main function for calculating BD shift.

    Parameters
    ----------
    args : dict
        See ``BD_shift`` subcommand
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
            msg = '  Comparing libraries: "{}", "{}"\n'
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
                            

def opt_parse(args=None):   
    docs = """
BD_shift: Determine the shift in BD based on KDE overlap

Usage:
  BD_shift [options] <kde1> <kde2>
  BD_shift -h | --help
  BD_shift --version

Options:
  <kde1>         KDE object
  <kde2>         KDE object 
  --start=<s>    Start of series for integration.
                 [default: 1.66]
  --end=<e>      End of series for integration.
                 [default: 1.85]
  --step=<x>     Step size of series for integration.
                 [default: 0.001]
  --np=<np>      Number of parallel processes.
                 [default: 1]
  --version      Show version.
  --debug        Debug mode (no parallel processes)
  -h --help      Show this screen.

Description:
  Determine the shift in BD value distribution between 2
  KDEs of BD values.
  The BD shift is calculated as 1 - KDE_intersection, 
  where KDE_intersection is calculated from 
  the integral of intersection densities at each point
  in the specified value series (start,end,step).
 
  Output
  ------
  A tab-delimited table of BD shift values is written to STDOUT.
    """
    if args is None:        
        args = docopt(docs, version='0.1')
    else:
        args = docopt(docs, version='0.1', argv=args)
    main(args)
          
