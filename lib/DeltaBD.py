"""Error distribution functions"""

# import
## batteries
import sys
from functools import partial
## 3rd party
import pandas as pd
import numpy as np
## application
from OTU_Table import OTU_table
import Utils



def evenly_spaced_BDs(BDs, n):
    """Interplate evenly spaced BDs
    BDs : iterable
       BDs for which to find min/max for setting BD range
    n : int
       Number of intervals (bins)
    """
    BDs = BDs.iloc[:,0].tolist()
    BD_min = min(BDs)
    BD_max = max(BDs)
    return np.linspace(BD_min, BD_max, n)


def center_of_mass(df, BD_intv):
    """Calculate the center of mass, which is the weighted mean BD,
    with weights as interpolated relative abundances.
    """
    rel_abund_intv = np.interp(BD_intv, df['BD_mid'], df['rel_abund'])
    if np.sum(rel_abund_intv) == 0:
        ri = np.random.randint(0, len(rel_abund_intv), 1)
        rel_abund_intv[ri] = 1e-20
    return np.average(BD_intv, axis=0, weights=rel_abund_intv)
    
def flattenHierarchicalCol(col, sep = '_'):
    """Flatten pandas multiindex columns
    col : pd.columns
    sep : str
    """
    if not type(col) is tuple:
        return col
    else:
        new_col = ''
        for leveli,level in enumerate(col):
            if not level == '':
                if not leveli == 0:
                    new_col += sep
                new_col += level
        return new_col


def deltaBD(Uargs):
    """Algorithm 
    # For each taxon:
    ## For each library:
    ### Linearly interpolate N evenly spaced relative abundance (RA) values 
    ### center_of_mass = weighted BD, where weights are interpolated RAs
    ## calc mean center_of_mass for treatments & controls
    ## deltaBD = mean_center_mass_treat - mean_center_mass_control

    Parameters
    ----------
    Uargs : dict
        See deltaBD.py
    """
    # loading OTU tables
    sys.stderr.write('Loading OTU tables...\n')    
    df_otu = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')
    exp_design = Utils.load_exp_design(Uargs['<exp_design>'])

    # setting evenly-spaced BDs for interpolation
    sys.stderr.write('Calculating delta BD...\n')        
    BD_intv = evenly_spaced_BDs(df_otu.select(['BD_mid']), int(Uargs['-b']))
    
    # Linear interpolation
    func = lambda x : center_of_mass(x, BD_intv)
    df_cm = df_otu.apply_by_group(func, 
                                  val_index='CM', 
                                  groups=['taxon', 'library'], 
                                  inplace=False)                                                    

    # Adding exp_design to OTU_table
    df_cm['library'] = df_cm['library'].astype(str)
    exp_design['library'] = exp_design['library'].astype(str)
    df_cm = df_cm.merge(exp_design, how='inner', on='library')

    # Calculating stdev CM for control or treatment
    func = lambda x : np.std(x['CM'])
    groups = ['taxon','sample_type']
    stdev_CM = df_cm.groupby(groups).apply(func).reset_index()[0]

    # Calculating mean CM for control or treatment
    func = lambda x : np.mean(x['CM'])
    groups = ['taxon','sample_type']
    df_cm = df_cm.groupby(groups).apply(func).reset_index()
    df_cm.columns = groups + ['mean_CM']
    df_cm['stdev_CM'] = stdev_CM


    # making table wide for sample_type
    df_cm = df_cm.pivot(index='taxon', columns='sample_type').reset_index()
    df_cm.columns = df_cm.columns.map(flattenHierarchicalCol)
#    print df_cm

    # Calculating delta BD 
    func = lambda x : x['mean_CM_treatment'] - x['mean_CM_control']
    df_cm['delta_BD'] =  df_cm.groupby(['taxon']).apply(func).reset_index()[0]
        
    # return
    return df_cm 
