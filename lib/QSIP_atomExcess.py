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
    # x = row in OTU table dataframe
    rel_abunds = x['prop_abs_abund_frac']
    BD = x['BD_mid']
    W = np.average(BD, weights=rel_abunds)
    return W

def calc_wAverage_density(otu, groups):
    """Calculated the weighted average density (W_ij), where weights are defined
    by taxon proportional absolute abundance.
    
    Parameters
    ----------
    otu : OTU_Table.OTU_table object
    groups : list-like
        Groupings for applying taxon_abundance / total_abundance
    """
    return otu.apply_by_group(_calc_wAve_density, 'density', 
                              groups=groups, inplace=False)

def calc_mean_density(densities, exp_design):
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
    mean_densities.columns = ['taxon', 'control', 'treatment']
    return mean_densities

def calc_density_shift(mean_densities):
    """Calculate density shift (Z_i) between mean densities (W_LABi - W_LIGHTi).
    
    Parameters
    ----------
    mean_densities : pandas.DataFrame
        output of calc_mean_density
    """
    f = lambda x : x['treatment'] - x['control']
    mean_densities['BD_diff'] = mean_densities.apply(f, axis=1)

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
        M = -0.4987282 * M_light2GC(M_light) + 9.974564 + M_light
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



def qSIP_atomExcess(Uargs):
    """METHOD (Calculate % atom incorp)
    # NOTE: need to specify which libraries are control vs treatment
    ## comm or incorp file?
    # per taxon:
    ## Calc: total gradient copy number per taxon
    ## Calc: taxon_density = weighted average of density 
    ###      (weights=copy numbers by fraction)
    ## Calc: mean_taxon_density (summed across rep gradients)
    ## Calc: density_shift = (mean_density__treat - mean_density__control)
    ## Calc: atom_fraction_excess = see Hungate et al., AEM

    Parameters
    ----------
    Uargs : dict
        See qSIP_atomExcess.py
    """
    # loading tables
    ## experimental design
    exp_design = pd.read_csv(Uargs['<exp_design>'], sep='\t', header=None)
    exp_design.columns = ['library', 'sample_type']
    
    ## OTU table 
    otu = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')
    
    #-- calculating BD shift (Z) --#
    # getting proportional absolute abundancs for each taxon
    groups = ('library', 'taxon',)
    calc_prop_abs_abund(otu, groups)
    
    # calculated weighted average BD (W_ij)        
    densities = calc_wAverage_density(otu, groups)
    
    # calculating means of W for control/treatment
    mean_densities = calc_mean_density(densities, exp_design)

    # calculating density shifts 
    calc_density_shift(mean_densities)

    #-- calculating atom excess --#
    ## control GC
    BD2GC_v = np.vectorize(BD2GC)
    f = lambda x : BD2GC_v(x['control'])
    mean_densities['control_GC'] = mean_densities.apply(f, axis=1)

    ## control molecular weight    
    GC2M_light_v = np.vectorize(GC2M_light)
    f = lambda x : GC2M_light_v(x['control_GC'])
    mean_densities['control_MW'] = mean_densities.apply(f, axis=1)
        
    ## max heavy molecular weight
    M_light2heavyMax_v = partial(M_light2heavyMax, isotope=Uargs['-i'])
    M_light2heavyMax_v = np.vectorize(M_light2heavyMax_v)
    f = lambda x : M_light2heavyMax_v(x['control_MW'])
    mean_densities['treatment_max_MW'] = mean_densities.apply(f, axis=1)
    
    ## molecular weight of labeled DNA
    calc_M_lab_v = np.vectorize(calc_M_lab)
    f = lambda x : calc_M_lab_v(x['BD_diff'], x['control'], x['control_MW'])
    mean_densities['treatment_MW'] = mean_densities.apply(f, axis=1)
        
    ## % atom excess
    # M_lab, M_light, M_heavyMax, isotope='13C'):
    calc_atomFracExcess_v = partial(calc_atomFracExcess, isotope=Uargs['-i'])
    calc_atomFracExcess_v = np.vectorize(calc_atomFracExcess_v)
    f = lambda x : calc_atomFracExcess_v(x['treatment_MW'], 
                                         x['control_MW'],
                                         x['treatment_max_MW'])
    mean_densities['atom_fraction_excess'] = mean_densities.apply(f, axis=1)
    
    return mean_densities
