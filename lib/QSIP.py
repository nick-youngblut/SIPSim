"""Error distribution functions"""

# import
## batteries
import sys
from functools import partial
## 3rd party
import numpy as np
import pandas as pd
## application
from OTU_Table import OTU_table


def neg_binom_err(m, r, negs=False):
    """Adding negative binomial distribuiton error, where variance
    scales more with the mean than a poisson distribution if (r < inf).

    Parameters
    ----------
    m : float
        Mean value
    r : float
        Negative binomial dispersion parameter
    negs : bool
        Negative values allowed? (otherwise 0)
    """
    sigma = np.sqrt(m + m**2 / r)
    x =  np.random.normal(m, sigma)
    if negs==False and x < 0:
        x = 0
    return x


def ave_neg_binom_err(m, r, negs=False, tech_reps=1):
    """Wrapper for neg_binom_err(). which will find the mean
    of multiple replicate error calculations (qPCR technical reps).

    Parameters
    ----------
    tech_reps : int
        Number of qPCR technical replicates.
    other_parameters : NA
        See neg_binom_err
    """
    errs = [neg_binom_err(m, r, negs) for i in xrange(tech_reps)]
    return np.mean(errs)

    
def write_qPCR(df, outFile):
    """Write out qPCR value table.
    
    Parameters
    ----------
    df : pandas dataframe
    outFile : string
        output file name
    """
    with open(outFile, 'wb') as outFH:
        df.to_csv(outFH, sep='\t', index=False)
    sys.stderr.write('qPCR value file written: {}\n'.format(outFile))


def qSIP(Uargs):
    """METHOD (apply qPCR data)
    # get total 'true' sample abundances from OTU table
    # calculate qPCR values
    ## get value(s) from neg-binom distribution
    ### `from Error_Dist import error_dist`
    #### alpha = ?
    ## [aside] write qPCR value table
    ## get mean of replicates 
    ### used to extrapolate abundances
    # transform rel-abundances to proportional abundances of total copy number
    # write table of transformed values 

    Parameters
    ----------
    Uargs : dict
        See qSIP.py

    """
    # loading OTU table (s)
    otu_abs = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')
    otu_rel = OTU_table.from_csv(Uargs['<OTU_subsample_table>'], sep='\t')


    # getting error from total OTU counts
    r = float(Uargs['-r'])
    reps = int(Uargs['--reps'])
    f = lambda x : ave_neg_binom_err(sum(x), r=r, tech_reps=reps)
    otu_abs.transform_by_group(f, ['count'], ['total_qPCR_copies'])


    # adding 'qPCR' values to subsampled (relative abundance) OTU table
    cols = ['library', 'fraction', 'taxon', 'total_qPCR_copies']
    otu_rel.merge(otu_abs.select(cols), how='inner', 
                  on=['library','fraction','taxon'])


    # Determining the proportional absolute abundances (y)
    ## y = total_16S_copies * taxon_rel_abundance
    groups = ['library', 'fraction', 'taxon']
    def _prop_abund(x):
        y = x['rel_abund'] * x['total_qPCR_copies']
        return  np.round(y, 0).astype(int)
    otu_rel.apply_by_group(_prop_abund, 'prop_abs_abund', groups)

    # return
    return otu_rel

