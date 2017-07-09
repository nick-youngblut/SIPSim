"""Error distribution functions"""

# import
## batteries
import sys
from functools import partial
## 3rd party
import pandas as pd
import sympy as sy
import numpy as np
## application
from SIPSim.OTU_Table import OTU_table

    
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


def calc_qPCR(m, expr, negs=False, tech_reps=1):
    """Add simulated error to a value (error represented by expr)
    to simulate qPCR copy number value(s).
    The mean of all technical replicates (tech_reps) will be returned.    

    Parameters
    ----------
    m : float        
    expr : lambdified sympy expression
        Expression for relating variance to the mean (var ~ mean)
    negs : bool
        Are negative values allowed?
    tech_reps : int
        Number of qPCR technical replicates.
    """
    
    var = expr(m)
    sigma = np.sqrt(int(var))
    x = np.random.normal(loc=m, scale=sigma, size=tech_reps)
    x = np.mean(x)
    if negs==False and x < 0:
        x = 0
    return x


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
    # loading OTU tables
    sys.stderr.write('Loading OTU tables...\n')    
    otu_abs = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')
    otu_rel = OTU_table.from_csv(Uargs['<OTU_subsample_table>'], sep='\t')

    # setting error function
    sys.stderr.write('Creating qPCR values...\n')    
    ## expression
    expr = sy.sympify(Uargs['-f'])
    x = sy.symbols('x')
    expr = sy.lambdify(x, expr, 'numpy')
    ## transforming values (first: summing all counts for each taxon w/ np.sum)
    tech_reps = int(Uargs['--reps'])
    f = lambda x : calc_qPCR(np.sum(x), expr=expr, tech_reps=tech_reps)
    otu_abs.transform_by_group(f, ['count'], ['total_qPCR_copies'])
    
    # adding 'qPCR' values to subsampled (relative abundance) OTU table
    sys.stderr.write('Calculating proportional absolute abundances...\n')    
    cols = ['library', 'fraction', 'taxon', 'total_qPCR_copies']
    otu_rel.merge(otu_abs.select(cols), how='inner', 
                  on=['library','fraction','taxon'])

    # Determining the proportional absolute abundances (y)
    ## y = total_16S_copies * taxon_rel_abundance
    otu_rel.df['prop_abs_abund'] = otu_rel.df['rel_abund'] * otu_rel.df['total_qPCR_copies']
    otu_rel.df['prop_abs_abund'] = otu_rel.df['prop_abs_abund'].astype(int)
    otu_rel.sort_values(by=['library', 'taxon', 'BD_mid'], inplace=True)

    # return
    return otu_rel

