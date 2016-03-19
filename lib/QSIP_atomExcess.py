"""Error distribution functions"""

# import
## batteries
import sys
from functools import partial
## 3rd party
import numpy as np
import pandas as pd
from pathos.multiprocessing import ProcessingPool
## application
from OTU_Table import OTU_table


def _prop_abund(x):
    """Calculate proportional absolute abundances.
    x = row in OTU dataframe
    """
    try:
        x = x['prop_abs_abund']
    except KeyError:
        msg = '"prop_abs_abund" column not found!' + \
              ' Check out qPCR subcommand.'
        sys.exit(msg)
    return (x / float(sum(x))).fillna(0)

def calc_prop_abs_abund(otu, groups):
    """Calculate the proporitonal absolute abundance (see Hungate et al., 2015)
    for each group specified by `groups`. 
    Prop_abs_abund = taxon copy number per sample as a fraction of total 
    copies for the taxon.
    
    Parameters
    ----------
    otu : OTU_Table.OTU_table object
    groups : list-like
        Groupings for applying taxon_abundance / total_abundance
    """
    otu.apply_by_group(_prop_abund, 'prop_abs_abund_frac', groups=groups)

def _calc_wAve_density(x):
    """Calculate weighted avareage density.
    x : row in OTU table dataframe
    """
    rel_abunds = x['prop_abs_abund_frac']
    BD = x['BD_mid']
    if np.sum(rel_abunds) <= 0:
        return np.nan
    else:
        return np.average(BD, weights=rel_abunds)

def calc_wAverage_density(otu, groups):
    """Calculated the weighted average density (W_ij), where weights are 
    defined by taxon proportional absolute abundance.
    
    Parameters
    ----------
    otu : OTU_Table.OTU_table object
    groups : list-like
        Groupings for applying taxon_abundance / total_abundance
    """
    return otu.apply_by_group(_calc_wAve_density, 'density', 
                              groups=groups, inplace=False)

def calc_mean_density(densities, exp_design=None):
    """Calculate mean densities (W_LIGHTi & W_LABi) for control/treatment 
    libraries.

    Parameters
    ----------
    densities : pandas.DataFrame
        output of calc_wAverage_density
    exp_design : pandas.DataFrame
        2-column dataframe (library, sample_type), 
          where sample_type = 'control' or 'treatment'
    """
    # calculating means of W for control/treatment
    ## join densities with exp_design
    if exp_design is not None:
        densities['library'] = densities['library'].astype(str)
        exp_design['library'] = exp_design['library'].astype(str)
        densities = pd.merge(densities, exp_design, how='inner', on='library')
    ## groupby, then apply (mean)
    groups = ['taxon', 'sample_type']
    f = lambda x : np.mean(x['density'])
    mean_densities = densities.groupby(groups).apply(f).reset_index()
    ncol = len(mean_densities.columns)
    mean_densities.columns = mean_densities.columns[:ncol-1].tolist() + \
                             ['mean_density']

    # wide to long table format
    mean_densities = mean_densities.pivot('taxon', 
                                          'sample_type', 
                                          'mean_density').reset_index()
    if mean_densities.shape[1] != 3:
        print mean_densities
        assert mean_densities.shape[1] == 3
    mean_densities.columns = ['taxon', 'control', 'treatment']
    return mean_densities

def calc_density_shift(df):
    """Calculate density shift (Z_i) between mean densities (W_LABi - W_LIGHTi).
    
    Parameters
    ----------
    df : pandas.DataFrame
        output of calc_mean_density
    """
    f = lambda x : x['treatment'] - x['control']
    df['BD_diff'] = df.apply(f, axis=1)

def BD2GC(BD):
    """buoyant density to G+C (fraction). 
    Ref: Birnie and Rickwood, 1978
    """
    return (BD - 1.66) / 0.098

def BD2GC_hungate(BD):
    """buoyant density to G+C (fraction). 
    Ref: Hungate et al., 2015. AEM
    """    
    return (1/0.083506) * (BD - 1.646057)

def GC2M_light(GC):
    """G+C to M_light (molecular weight of 'light' DNA).
    Ref: Hungate et al., 2015. AEM
    """
    return 0.496 * GC + 307.691

def M_light2GC(M_light):
    """The inverse of GC2M_light.
    """
    return (M_light - 307.691) / 0.496
    
def M_light2heavyMax(M_light, isotope='13C'):
    """G+C to theoretical molecular weight of fully-labeled DNA (M_HEAVYMAXi).
    Ref: Hungate et al., 2015. AEM.
    
    Parameters
    ----------
    M_light : float
        Molecular weight of DNA in control gradients (light DNA)
    isotope : str
        13C or 18O isotope?
    """
    if isotope == '13C':        
        G = M_light2GC(M_light)
        M = -0.4987282 * G + 9.974564 + M_light
    elif isotope == '18O':
        M = 12.07747 + M_light
    else:
        msg = '"{}" Isotope not supported; only: "13C" or "18O"'
        raise TypeError, msg.format(isotope)
    return M
    
def calc_M_lab(Z, W_light, M_light):        
    """Calculate the molecular weight of DNA in labeled treatments.
    Ref: Hungate et al., 2015. AEM.
    
    Parameters
    ----------
    Z : float
        Difference in mean densities between labeled and control gradients
    W_light : float
        Mean density of all control gradients
    M_light : float
        Molecular weight of DNA in control gradients
    """
    return (Z / W_light + 1) * M_light 

def calc_atomFracExcess(M_lab, M_light, M_heavyMax, isotope='13C'):
    """Calculate atom fraction excess from qSIP data.
    
    Parameters
    ----------
    M_lab : float
        Molecular weight of DNA in labeled treatment gradients
    M_light : float
        Molecular weight of DNA in control treatment gradients       
    M_heavyMax : float
        Theoretical molecular weight of fully-labeled DNA
    isotope : str
        13C or 18O isotope?
    """
    if isotope == '13C':
        a = 0.01111233
    elif isotope == '18O':
        a = 0.002000429
    else:
        msg = '"{}" Isotope not supported; only: "13C" or "18O"'
        raise TypeError, msg.format(isotope)

    x = M_lab - M_light
    y = M_heavyMax - M_light
    return x / y * (1 - a)


def atomFracExcess(mean_densities, isotope='13C'):
    """Calculate atom fraction excess.

    Parameters
    ----------
    mean_densities : pandas.DataFrame
        Dataframe of weighted mean densities averaged across gradients.
    isotope : string
        Which isotope to calculate?
    """
    ## control GC
    BD2GC_v = np.vectorize(BD2GC)
    f = lambda x : BD2GC_v(x['control'])
    mean_densities['control_GC'] = mean_densities.apply(f, axis=1)

    ## control molecular weight    
    GC2M_light_v = np.vectorize(GC2M_light)
    f = lambda x : GC2M_light_v(x['control_GC'])
    mean_densities['control_MW'] = mean_densities.apply(f, axis=1)
        
    ## max heavy molecular weight
    M_light2heavyMax_v = partial(M_light2heavyMax, isotope=isotope)
    M_light2heavyMax_v = np.vectorize(M_light2heavyMax_v)
    f = lambda x : M_light2heavyMax_v(x['control_MW'])
    mean_densities['treatment_max_MW'] = mean_densities.apply(f, axis=1)
    
    ## molecular weight of labeled DNA
    calc_M_lab_v = np.vectorize(calc_M_lab)
    f = lambda x : calc_M_lab_v(x['BD_diff'], x['control'], x['control_MW'])
    mean_densities['treatment_MW'] = mean_densities.apply(f, axis=1)
        
    ## % atom excess
    calc_atomFracExcess_v = partial(calc_atomFracExcess, isotope=isotope)
    calc_atomFracExcess_v = np.vectorize(calc_atomFracExcess_v)
    f = lambda x : calc_atomFracExcess_v(x['treatment_MW'], 
                                         x['control_MW'],
                                         x['treatment_max_MW'])
    mean_densities['atom_fraction_excess'] = mean_densities.apply(f, axis=1)


def subsample_densities(df, sample_type):
    """Subsampling density values for bootstrapping.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame that contains 'sample_type'  (values: control or treatment)
    sample_type : string
        Which sample type (control|treatment) to sub-sample for.
    """
    size = df.loc[df['sample_type'] == sample_type].shape[0]
    subsamps = np.random.choice(df.index, size, replace=True)
    x = df.loc[subsamps,:].reset_index()
    x['sample_type'] = sample_type
    return x


def _bootstrap(df, isotope):
    # subsample with replacement
    ## control
    idx = df['sample_type'] == 'control'
    df_cont = subsample_densities(df, 'control')
    ## treatment
    idx = df['sample_type'] == 'treatment'
    df_treat = subsample_densities(df, 'treatment')
    ## cat
    df_i = pd.concat([df_cont, df_treat])
    
    # calculate weighted-mean density (wmd)
    df_i_wmd = calc_mean_density(df_i)
    
    # calculating density shifts 
    calc_density_shift(df_i_wmd)
    
    # calculating atom fraction excess 
    atomFracExcess(df_i_wmd, isotope=isotope)
    
    # return
    assert df_i_wmd.shape[0] == 1
    return df_i_wmd.loc[0,:]


def _bootstrap_CI(df, n=1000, a=0.1, isotope='13C', pool=None):
    # status
    taxon = df['taxon'].unique()[0]
    msg = 'Bootstrap CI (n={}); processing taxon: {}\n'
    sys.stderr.write(msg.format(n, taxon))

    # bootstrapping
    if pool is None:        
        boot_res = [_bootstrap(df, isotope) for i in xrange(n)]
    else:
        dfs = [df for i in xrange(n)]
        _bootstrap_p = partial(_bootstrap, isotope=isotope)
        boot_res = pool.map(_bootstrap_p, dfs)
    boot_res = pd.concat(boot_res, axis=1) 
    boot_res.columns = range(boot_res.shape[1])
    boot_res = boot_res.transpose()
        
    # calculating deltas [true_value - bootstrap_value] 
    true_A = df.reset_index().loc[0,'atom_fraction_excess']
    f = lambda x : true_A - x['atom_fraction_excess']
    A_delta = boot_res.apply(f, axis=1).tolist()
    
    # calculating CI
    perc_low = np.percentile(A_delta, a / 2 * 100)
    perc_high = np.percentile(A_delta, (1 - a / 2) * 100)
    if perc_low > true_A:
        perc_low = true_A
    if perc_high < true_A:
        perc_high = true_A
    CI_low = true_A - np.abs(perc_low)
    CI_high = true_A + np.abs(perc_high)
    CI_low = np.round(CI_low, 6)
    CI_high = np.round(CI_high, 6)
    msg = 'WARNING: CI_low ({}) is > CI_high ({})\n'
    if CI_low > CI_high:
        sys.stderr.write(msg.format(CI_low, CI_high))
    
    return pd.Series({'atom_CI_low' : CI_low, 'atom_CI_high' : CI_high})


def bootstrap_CI(densities, mean_densities, exp_design,
                 n=1000, a=0.1, isotope='13C', nodes=1):
    """Calculate qSIP bootstraped confidence intervals.
    Reference: Hungate et al., 2015.

    Parameters
    ----------
    densities : pandas.DataFrame
        Table of all weighted mean densities for each library-taxon.
    mean_densities : pandas.DataFrame
        Table of library-averaged densities for each taxon.
    exp_design : pandas.DataFrame
        Table with 'library' and 'sample_type' columns (control|treatment)
    n : int
        Number of bootstrap replicates
    a : float
        Alpha for confidence interval calculation
    isotope : str
        Which isotope to calculate atom fraction excess?
    nodes : int
        Number of parallel processes.    
    """    
    # multiprocessing
    if nodes is None:
        pool = None
    else:
        pool = ProcessingPool(nodes=nodes)

    # add: mean_densities
    cols = ['taxon', 'BD_diff', 'atom_fraction_excess']
    densities = densities.merge(mean_densities[cols], on=['taxon'], how='inner')
    # add: experimental design
    densities = densities.merge(exp_design, on=['library'], how='inner')

    # calculate CI
    f = lambda x : _bootstrap_CI(x, n=n, a=a, isotope=isotope, pool=pool)
    CIs = densities.groupby(['taxon']).apply(f).reset_index()
    cols = ['taxon', 'atom_CI_low', 'atom_CI_high']
    return CIs[cols]


def load_exp_design(inFile):
    exp_design = pd.read_csv(inFile, sep='\t', header=None)
    exp_design.columns = ['library', 'sample_type']
    # formatting
    f = lambda x : x.lower()
    exp_design['sample_type'] = exp_design['sample_type'].apply(f)
    # assert 
    x = exp_design['sample_type'].isin(['control','treatment'])
    msg = 'Only ("control" or "treatment" allowed in 2nd-column' + \
          'of <exp_design> table'
    assert all(x) == True, msg
    # return
    return exp_design


def qSIP_atomExcess(Uargs):
    """Main function for calculating atom fraction excess (and density shift)
    for qSIP data (OTU table where relative abundances were multipled by qPCR
    values (total gene copies) to get proportinoal absolute abundances).

    Parameters
    ----------
    Uargs : dict
        See qSIP_atomExcess.py
    """
    # loading tables
    sys.stderr.write('Loading files...\n')
    ## experimental design
    exp_design = load_exp_design(Uargs['<exp_design>'])
        
    ## OTU table 
    otu = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')
    
    #-- calculating BD shift (Z) --#
    sys.stderr.write('Calculating density shifts (Z)...\n')
    # getting proportional absolute abundancs for each taxon
    groups = ('library', 'taxon',)
    calc_prop_abs_abund(otu, groups)
    
    # calculated weighted average BD (W_ij)        
    densities = calc_wAverage_density(otu, groups)
    
    # calculating means of W for control/treatment
    mean_densities = calc_mean_density(densities, exp_design)

    # calculating density shifts 
    calc_density_shift(mean_densities)

    #-- calculating atom fraction excess --#
    sys.stderr.write('Calculating atom fraction excess (A)...\n')
    atomFracExcess(mean_densities, isotope=Uargs['-i'])

    #-- calculating CIs --#
    sys.stderr.write('Calculating bootstrap CIs...\n')
    if Uargs['--debug']:
        taxa = mean_densities['taxon'].tolist()
        densities = densities.loc[densities['taxon'].isin(taxa[:20])]
        Uargs['--np'] = None
    else:
        Uargs['--np'] = int(Uargs['--np'])
    CIs = bootstrap_CI(densities, mean_densities, exp_design, 
                       n=int(Uargs['-n']), 
                       a=float(Uargs['-a']), 
                       nodes=Uargs['--np'],
                       isotope=Uargs['-i'])
    mean_densities = mean_densities.merge(CIs, on=['taxon'], how='inner')
    
    # return
    return mean_densities
