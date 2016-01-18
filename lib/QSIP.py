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

    # getting total absolute abunds for each sample
    f = lambda x : sum(x['count'])
    total_cnt = otu_abs.apply_by_group(f, 'total_count', 
                                       inplace=False).reset_index()
    
    # drawing error from OTU counts
    f = partial(neg_binom_err, r=float(Uargs['-r']))
    f = np.vectorize(f)
    total_cnt['total_count_qPCR'] =  f(total_cnt['total_count'])

    # writing out total count table
    write_qPCR(total_cnt, Uargs['-o'])

    # Transforming relative abunds to abs abunds for each taxon (in each sample)
    ## determining the proportional absolute abundances (y)
    ## y = total_gene_copies * taxon_rel_abundance
    otu_rel.merge(total_cnt, how='inner', on=['library','fraction'])
    groups = ['library', 'fraction', 'taxon']
    def _prop_abund(x):
        y = x['rel_abund'] * x['total_count_qPCR']
        return  np.round(y, 0).astype(int)
    otu_rel.apply_by_group(_prop_abund, 'prop_abs_abund', groups)

    # deleting 'temporary' calculation columns
    otu_rel.drop('index',1)
    otu_rel.drop('total_count',1)
    otu_rel.drop('total_count_qPCR',1)

    return otu_rel

